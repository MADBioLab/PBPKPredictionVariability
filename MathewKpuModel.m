function Params=MathewKpuModel(Params)
    %% read in input parameters
    if strcmp(Params.species,'hum')
        physparams=importdata('tissueCompositionParamsHuman.csv');
    elseif strcmp(Params.species,'rat')
        physparams=importdata('tissueCompositionParamsRat.csv');
    end

    tissues=physparams.textdata(2:end,1);
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

    %normalize tissue fractions
    Total=fNL+fNP+fEW+fIW+APfrac;
    fNL=fNL./Total;
    fNP=fNP./Total;
    fEW=fEW./Total;
    fIW=fIW./Total;
    APfrac=APfrac./Total;

    %% Impose lower bound on Kpu
    Params.fu=max([Params.fu,1e-4]);
    
    %% physiological pH values
    pHIW = 7.0;        %intracellular water
    pHp = 7.4;         %plasma and extracellular water
        
    %% Perform ionization calculations for plasma/extracellular and intracellular pH
    fI = SimCypIon(Params.pKa_a,Params.pKa_b,pHp); 
    fI_IW = SimCypIon(Params.pKa_a,Params.pKa_b,pHIW);
    Params.chargeType = fI.chargeType;
    Params.netCharge = fI.netCharge;

    %% logP calculations (including logD and logPvow conversions)    
    % infer logP from logD74 and ionization microstates
    Params.logP=calcLogP(fI,Params.logD74);
    P=10^Params.logP;

    % for adipose: estimate vegetable oil-water partitioning from octanol
    % water partitioning
    logPvow=1.115*Params.logP-1.35;
    Pvow = 10^logPvow;    
    
    %% Define function to calculate total neutral lipid and neutral phospholipid association
    KnplFunc = @(P) (0.3*P+0.7); 
    lipidBindFunc=@(fNL,fNP,P) fNL*P + fNP*KnplFunc(P);
    
    %% Approximate protein binding affinity 
    Ka = 1/(Params.fu+eps)-fEW(end);    

    %% Adjust fraction unbound to account for plasma lipids
    temp = Ka*((1-fI.fN)*Alb(end)+fI.fN*lipoprot(end)) ...
         + fI.fN * lipidBindFunc(fNL(end),fNP(end),P);
    Params.fu_adj=1-temp/(temp+fEW(end));
    
    %% Calculate tissue partitioning for each organ and red blood cells 
    Params.Kpu=zeros(Ntissues,1);    
    Params.IW2EW=zeros(Ntissues,1);

    for j=[1:Ntissues-2,Ntissues]
        P=10^Params.logP;
        %% use vegetable oil-water partitioning for Adipose
        flag = strcmp(tissues{j}, 'Adipose'); 
        P = flag * Pvow + ~flag*P;        
        Kapl = KaplFunc(fI_IW.fN,0,fI_IW.fA,fI_IW.fC,KnplFunc(P));

        %% total tissue partitioning
        Params.IW2EW(j) = fI.fN/(fI_IW.fN+eps);         
        Params.Kpu(j) = fEW(j)+Ka*((1-fI.fN)*Alb(j) + fI.fN * lipoprot(j)) ...
                        + Params.IW2EW(j) * (fIW(j)+fI_IW.fC*Kapl*APfrac(j)) ...
                        +  fI.fN * lipidBindFunc(fNL(j),fNP(j),P); 
    
        Params.(strcat('Kpu_',(tissueNames{j}))) = Params.Kpu(j);
      
    end
    Params.CeCp = Params.Kpu_RBCs * Params.fu_adj;
    Params.B2P = Params.H * Params.CeCp+1 - Params.H;

end

%%

function f=calcLogP(fI,logD)
    dlogPmc = 3;            
    dlogPma = 4;            
    dlogPdc = 5;             
    dlogPda = 5;            
    dlogPzw = 2.5;           
    dlogPzw2 = 5;            
    IPP=1;
      
    f = logD - log10(fI.fN + IPP * (fI.fb / 10^dlogPmc + fI.fbb / 10^dlogPdc ...
        + fI.fa / 10^dlogPma + fI.fab / 10^dlogPzw + fI.fabb / 10^dlogPdc ...
        + fI.faa / 10^dlogPda + fI.faab / 10^dlogPda + fI.faabb / 10^dlogPzw2)) / log10(10);    
       
end

%%

function fI = SimCypIon(pKa_a,pKa_b,pH)
    %% classify whether site should be considered ionizable 
    pKAcidLim=9.4; %acidic site if pKa_a<9.4
    pKBaseLim=5.4; %basic site if pKa_b>5.4
    
    %% classify ion class of molecule
    if pKa_a<=pKAcidLim
        fI.chargetype='acid';
        if pKa_b>pKBaseLim
            fI.chargetype='zwit';
        end
    elseif pKa_b>pKBaseLim
        fI.chargetype='base';
    else
        fI.chargetype='neut';
    end
    
    KaA1 = 10 ^ (-pKa_a(1));
    KaA2 = 10 ^ (-pKa_a(2));
    KaB1 = 10 ^ (-pKa_b(1));
    KaB2 = 10 ^ (-pKa_b(2));
    H = 10 ^ (-pH);
    denom = ((H + KaA1) * (H + KaA2) * (H + KaB1) * (H + KaB2));
    %neutral
    fI.fN = H ^ 2 * KaB1 * KaB2 / denom;%                      'unionized (neutral)
    %cations
    fI.fb = H ^ 3 * (KaB1 + KaB2) / denom;%                    'mono-cation (+1)
    fI.fbb = H ^ 4 / denom;%                                   'di-cation (+2)
    fI.fabb = H ^ 3 * (KaA1 + KaA2) / denom;%                  'mono-anion & di-cation (+1)
    %anions
    fI.fa = H * KaB1 * KaB2 * (KaA1 + KaA2) / denom;%          'mono-anion (-1)
    fI.faa = KaA1 * KaA2 * KaB1 * KaB2 / denom;%               'di-anion (-2)
    fI.faab = H * KaA1 * KaA2 * (KaB1 + KaB2) / denom;%        'di-anion & mono-cation (-1)
    %zwitterions
    fI.fab = H ^ 2 * (KaA1 + KaA2) * (KaB1 + KaB2) / denom;%   'mono-anion & mono-cation (neutral)
    fI.faabb = H ^ 2 * KaA1 * KaA2 / denom;%                   'di-anion & di-cation (neutral)
    
    fI.fA = fI.fa+fI.faa;
    fI.fC = fI.fb+fI.fbb+fI.fabb + fI.faab+fI.fab+fI.faabb; %all carriers of a protonated site
    fI.netCharge = (fI.fb+fI.fabb) + 2*fI.fbb - (fI.fa+fI.faab) - 2*fI.faa;
    
end
