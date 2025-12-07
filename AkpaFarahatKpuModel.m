function Params=AkpaFarahatKpuModel(Params)

    %% read in tissue composition parameters
    if strcmp(Params.species,'hum')
        physparams=importdata('tissueCompositionParamsHuman.csv');
    elseif strcmp(Params.species,'rat')
        physparams=importdata('tissueCompositionParamsRat.csv');
    end

    physparams=physparams.data;
    tissueNames=[Params.tissueNames;'RBCs'];
    Ntissues=length(tissueNames); 

    %fractional tissue volumes occupied by:
    fNL=physparams(:,1);        %neutral lipids
    fNP=physparams(:,2);        %neutral phospholipids
    fEW=physparams(:,3);        %extracellular water
    fIW=physparams(:,4);        %intracellular water    
    lipoprot=physparams(:,5);   %lipoproteins
    Alb=physparams(:,6);        %albumin
    APfrac=physparams(:,7);     %acidic phospholipids
    ICprot=physparams(:,8);     %intracellular proteins

    %normalize tissue fractions
    Total=fNL+fNP+fEW+fIW+APfrac+ICprot;
    fNL=fNL./Total;
    fNP=fNP./Total;
    fEW=fEW./Total;
    fIW=fIW./Total;
    APfrac=APfrac./Total;
    ICprot=ICprot./Total;

    %tissue-specific intracellular pH
    pH=physparams(:,9);
    
    %% Impose lower bound of 1e-4 on fu
    Params.fu=max([Params.fu,1e-4]);
    
    %% Perform ionization calculations for plasma/extracellular pH
    [Params.chargeType,Params.netCharge,I_PL]=...
                          ionization(Params.pKa_a,Params.pKa_b,pH(end));   

    %% logP calculations (including logD and logPvow conversions)    
    D74 = 10^Params.logD74; 
    % infer logP from logD74 and ionization microstates
    Params.logP=calcLogP(I_PL.fn,Params.logD74);
    % for adipose: estimate vegetable oil-water partitioning from octanol
    % water partitioning
    logPvow=1.115*Params.logP-1.35;
    Pvow = 10^logPvow;    
    
    %% Define function to estimate neutral phospholipid association
    %membrane affinity correlation of Yun and Edginton
    KnplFunc = @(D) 10^(1.294+0.304*log10(D)); 

    %% Define function to calculate total neutral lipid and neutral phospholipid association
    lipidBindFunc=@(fNL,fNP,Plip,Pplip) fNL*Plip + fNP*KnplFunc(Pplip);

    %% Define function to calculate acidic phospholipid partition coefficient
    % fn, fz, fa, fc are neutral, zwitterionic, anionic, and cationic
    % fractions, respectively
    KaplFunc = @(fn,fz,fa,fc,Knpl) Knpl*(fn + 0.05 * fa + 20*(fc+fz));

    %% Define function to estimate intracellular protein binding
    % D is the intracellular distribution coefficient
    KprotFunc=@(D) 0.163 + 0.0221*KnplFunc(D);

    %% 'Typical' protein binding coefficient, estimated from fu
    KaPRp = fEW(end)*(1/(Params.fu+eps)-1); %Kapr*PRplasma per unit volume 

    %% Adjust fraction unbound to account for plasma lipids
    % Calculate Kapl at plasma pH
    Kapl = KaplFunc(I_PL.fn+eps,I_PL.fz,I_PL.fa,I_PL.fc,KnplFunc(10^Params.logP));
    Params.fu_adj = fEW(end)/(fEW(end) + KaPRp ...
           + lipidBindFunc(fNL(end),fNP(end),D74,D74) ...
           + Kapl*APfrac(end));  
   
    %% Calculate tissue partitioning for each organ and red blood cells 
    Params.Kpu=zeros(Ntissues,1);    
    Params.IW2EW=zeros(Ntissues,1); %ratio of intracellular to extracellular unbound concentrations

    for j=[1:Ntissues-2,Ntissues] 
        %% Ionization calculations for intracellular pH
        [~,~,I]=ionization(Params.pKa_a,Params.pKa_b,pH(j));

        %% use vegetable oil-water partitioning for Adipose
        flag = strcmp(tissueNames{j}, 'Adipose'); 
        P = flag * Pvow + ~flag * 10^Params.logP;
        
        %% intracellular distribution coefficient (i.e., @intracellular pH)
        % estimate distribution coefficient at intracellular pH
        DIW=10^calcLogD(I.fn,log10(P)); 
              
        %% unbound in intracellular water
        term1 = fIW(j);
       
        %% association with intracellular neutral lipids, neutral phospholipids, and proteins
        term2 = lipidBindFunc(fNL(j),fNP(j),DIW,DIW) ...
              + KprotFunc(10^Params.logP) * ICprot(j); 
        
        %% association w/ proteins in extracellular water 
        term3 = (I_PL.fn * lipoprot(j) + (I_PL.fa+I_PL.fz) * Alb(j)) * KaPRp;        
        
        %% association with acidic phospholipids  
        Kapl = KaplFunc(I.fn,I.fz,I.fa,I.fc,KnplFunc(10^Params.logP));
        term4 = Kapl*APfrac(j);

        %% total tissue partitioning
        Params.IW2EW(j) = I_PL.fn/(I.fn+eps);     
        Params.Kpu(j) = (fEW(j) + term3) + Params.IW2EW(j) * (term1 + term2 + term4);

        Params.(strcat('Kpu_',(tissueNames{j}))) = Params.Kpu(j);
        
    end

    %% calculate B2P
    %ratio of concentration in RBCs to total concentration in plasma
    Params.KRBC = Params.Kpu_RBCs*Params.fu_adj; 
    Params.B2P = Params.H * KRBC+1-Params.H;

end

%%
function f=calcLogD(fn,logP)
    %use ionization state and logP to estimate logD
    ratio=fn+1e-3*(1-fn);
    f=log10(ratio*10^logP);
end

%%
function f=calcLogP(fn,logD)
    %use ionization state and logD to estimate logP
    ratio=fn+1e-3*(1-fn);
    f=log10(10^logD/ratio);
end
%%

function [chargetype,netCharge,I]=ionization(pKa_a,pKa_b,pH)
    %% classify whether site should be considered ionizable 
    pKAcidLim=9.4; %acidic site if pKa_a<9.4
    pKBaseLim=5.4; %basic site if pKa_b>5.4
    
    %% classify ion class of molecule
    if pKa_a<=pKAcidLim
        chargetype='acid';
        if pKa_b>pKBaseLim
            chargetype='zwit';
        end
    elseif pKa_b>pKBaseLim
        chargetype='base';
    else
        chargetype='neut';
    end
    
    %% determine ionization probability for each ionizable site
    ionzn=zeros(1,2);
    pKa=[pKa_a,pKa_b];
    ionizationType={'acid','base'};
    for i=1:length(pKa)
        ionzn(i)=HendersonHasselbalch(ionizationType{i},pH,pKa(i))/100; % [%]
    end
    
    %% consider all possible combinations of charged sites and calculate
    %% weighted molecular charge contribution (weighted by probability of each
    %% molecular ionization state)
    Pcharged=ionzn;
    Puncharged=100-ionzn;
    % -1 charge
    fa_combinations = [1 0];
    % +1 charge
    fc_combinations = [0 1];
    % neutral
    fn_combinations = [0 0];
    % zwitterions
    fz_combinations = [1 1];
    % all molecules with a protonated site 
    fwithPos_combinations = [[0 1];
                             [1 1]];
    
    %calculate fraction of molecule population with +1, -1, or 0 charge
    fc=fc_combinations.*Pcharged ...
        + ~fc_combinations.*Puncharged;
    fc=sum(prod(fc,2));

    fa=fa_combinations.*Pcharged ...
        + ~fa_combinations.*Puncharged;
    fa=sum(prod(fa,2));

    fn=fn_combinations.*Pcharged ...
        + ~fn_combinations.*Puncharged;
    fn=sum(prod(fn,2));

    fz=fz_combinations.*Pcharged ...
        + ~fz_combinations.*Puncharged;
    fz=sum(prod(fz,2));
    
    %fraction of molecule population bearing a protonated site
    fwithPos=fwithPos_combinations.*Pcharged ...
          + ~fwithPos_combinations.*Puncharged;
    fwithPos=sum(prod(fwithPos,2));
    
    netCharge=fc-fa;
    
    I.fc=fc;
    I.fa=fa;
    I.fn=fn;
    I.fz=fz;
    I.fwithPos=fwithPos;

end

