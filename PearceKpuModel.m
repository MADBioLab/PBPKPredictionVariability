function Params=PearceKpuModel(Params)

    %% read in input parameters
    if strcmp(Params.species,'hum')
        physparams=importdata('tissueCompositionParamsPearceHuman.csv');
    elseif strcmp(Params.species,'rat')
        physparams=importdata('tissueCompositionParamsPearceRat.csv');
    end
    regressionParams=importdata('PearceRegressionParams.csv'); %INTERCEPT, SLOPE
    regressionParams=regressionParams.data;

    physparams=physparams.data;
    tissueNames=[Params.tissueNames;'RBCs'];
    Ntissues=length(tissueNames); 

    %fractional tissue volumes occupied by:
    Fcell=physparams(:,1);          %cell
    Fint=physparams(:,2);           %interstitial
    FwCell=physparams(:,3);         %cell water
    FlipidCell=physparams(:,4);     %cell total lipid
    FprotCell=physparams(:,5);      %cell proteins
    FNLCell=physparams(:,6);        %cell neutral lipids (as fraction of total lipid)
    FNPLCell=physparams(:,7);       %cell phospholipids
    FAPLCell=physparams(:,8);       %cell acidic phospholipids
    pH=physparams(:,11);            %intracellular pH 

    %% Impose lower bound on fu
    Params.fu=max([Params.fu,1e-4]);

    %normalize tissue fractions
    Fcell=Fcell./(Fcell+Fint);
    Fint=Fint./(Fcell+Fint);

    %% Perform fractional volume conversion
    FNLCell=FNLCell.*FlipidCell;        %cell neutral lipids
    FNPLCell=FNPLCell.*FlipidCell;      %cell phospholipids
    FAPLCell=FAPLCell.*FlipidCell;      %cell acidic phospholipids    
    FprotInt=FprotCell(end)*0.37;       %interstitial protein
    FwPl=1-FprotCell(end);              %plasma water
    FwInt=FwPl;                         %interstitial water

    %% Perform ionization calculations for plasma/extracellular and intracellular pH
    [Params.chargeType,Params.netCharge,I_PL]=...
                        ionization(Params.pKa_a,Params.pKa_b,pH(end));     

    %% logP calculations (including logD and logPvow conversions)    
    % infer logP from logD74 and microstates
    D74 = 10^Params.logD74; 
    Params.logP=calcLogP(I_PL.fn,Params.logD74);  
    
    %% Define function to estimate neutral phospholipid association
    KnplFunc = @(P) 10^(1.294+0.304*log10(P)); 
    %% Define function to calculate total neutral lipid and neutral phospholipid association
    lipidBindFunc=@(fNL,fNP,Plip,Pplip) fNL*Plip + fNP*KnplFunc(Pplip);
    %% Define function to calculate acidic phospholipid partition coefficient
    KaplFunc = @(fn,fz,fa,fc,Knpl) Knpl * (fn + fz + 0.05*fa + 20*fc);
    %% Define function to estimate intracellular protein binding
    KprotFunc=@(P) 0.163 + 0.0221*KnplFunc(P);

    %% Adjust fraction unbound to account for plasma lipids
    Params.fu_adj = 1/(FNLCell(end) * D74 + 1/(Params.fu+eps));       

    %% Calculate tissue partitioning for each organ and red blood cells 
    Kint = (FwInt + FprotInt/FprotCell(end) * (1/(Params.fu_adj+eps) - FwPl));
    Params.Kpu=zeros(Ntissues,1);    
    Params.IW2EW=zeros(Ntissues,1);

    for j=[1:Ntissues-2,Ntissues]
        %% Ionization calculations for intracellular pH
        [~,~,I]=ionization(Params.pKa_a,Params.pKa_b,pH(j));
        P=10^Params.logP;
        
        %% intracellular partition coefficient (@intracellular pH)
        PIW=10^calcLogD(I.fn,log10(P));
        
        %% calculate Kpu (Plasma unbound)     
        %% unbound in intracellular water
        term1 = FwCell(j);
       
        %% association with intracellular neutral lipids, neutral phospholipids, and proteins
        term2 = lipidBindFunc(FNLCell(j),FNPLCell(j),PIW,KnplFunc(10^Params.logP))...
                            + FprotCell(j) * KprotFunc(10^Params.logP);

        %% association with acidic phospholipids  
        Kapl = KaplFunc(I.fn,I.fz,I.fa,I.fc,KnplFunc(10^Params.logP));        
        term3 = FAPLCell(j) * Kapl;        
        
        %% cell partitioning
        Kcell = term1 + term2 + term3;

        %% total tissue partitioning
        KAPPAcell2pu = (I_PL.fn + I_PL.fz + 1e-3*(I_PL.fc+I_PL.fa)) ...
                     / (I.fn + I.fz + 1e-3*(I.fc+I.fa));        
        Params.IW2EW(j) = KAPPAcell2pu; 
        Params.Kpu(j) = Fint(j) * Kint + Params.IW2EW(j) * Fcell(j) * Kcell;

        %% apply regression correction  
        if j~=Ntissues
            Params.Kpu(j) = 10^(regressionParams(j,2) ...
                              * log10(Params.Kpu(j) * Params.fu_adj) ...
                              + regressionParams(j,1))...
                              /Params.fu_adj;
        end
        Params.(strcat('Kpu_',(tissueNames{j}))) = Params.Kpu(j);
      
    end
    
    %% calculate B2P
    Params.KRBC = Params.Kpu_RBCs * Params.fu_adj;
    Params.B2P = Params.H * Params.KRBC + 1-Params.H;

end

%%

function logD=calcLogD(fn,logP)
    %use ionization state and logP to estimate logD
    ratio=fn+1e-3*(1-fn);
    logD=log10(ratio*10^logP);
end

%%

function logP=calcLogP(fn,logD)
    %use ionization state and logD to estimate logP
    ratio=fn+1e-3*(1-fn);
    logP=log10(10^logD/ratio);
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
