function [ ODEs ] = PBPKodes(time, m, P, Y, bodyWeight)

    %% arguments
    % time - current simulation time [s]
    % m - current total mass of drug in each organ [mg]
    % P - parameter struct
    % Y - struct containing organ subcompartment concentrations for all organs
    % bodyWeight - total bodyweight [kg]

    %% Reference Units
    % concentrations -> moles/L , Q -> [L/s]
    
    %% drug administration
    % converts IV dose in mg/kg -> mole/s using body weight [kg] and duration [secs]
    % seconds. RMM is relative molar mass of drug
    source = (time <= P.dose_duration) * bodyWeight ...
            * P.dose / P.dose_duration / P.RMM / 1e3; 
         
    %% solve non-linear algebriac relationships for drug distribution in 
    %% organ sub-compartments
    Y = calcAlgebraic(m,P,Y);
    
    %% evaluate accumulation rates for each organ or blood compartment
    ODEs = zeros(length(m),1);
    
    %% Artery %%
    ODEs(1) = P.Q_Lung * Y.C_vasc_Lung ...
              - Y.C_vasc_Artery * (P.Q_Gut ...
              + P.Q_Spleen + P.Q_Heart + P.Q_Muscle...
              + P.Q_Adipose + P.Q_Kidney + P.Q_Skin...
              + P.Q_Brain + P.Q_Bone + P.Q_Liver_hep ...
              + P.Q_RestOfBody);% - metabolism;
              
    %% Vein %%
    ODEs(2) = source + P.Q_Heart * Y.C_vasc_Heart ...
              + P.Q_Muscle * Y.C_vasc_Muscle ...
              + P.Q_Adipose * Y.C_vasc_Adipose ...
              + P.Q_Kidney * Y.C_vasc_Kidney ...
              + P.Q_Skin * Y.C_vasc_Skin  ...
              + P.Q_Brain * Y.C_vasc_Brain ...
              + P.Q_Bone * Y.C_vasc_Bone ...
              + P.Q_Liver * Y.C_vasc_Liver ...
              + P.Q_RestOfBody * Y.C_vasc_RestOfBody ...
              -  P.Q_Vein * Y.C_vasc_Vein;

    %% Lung %%
    ODEs(3) = P.Q_Vein * Y.C_vasc_Vein -  P.Q_Lung * Y.C_vasc_Lung;
    %% Gut %%
    ODEs(4) = P.Q_Gut * Y.C_vasc_Artery -  P.Q_Gut * Y.C_vasc_Gut;
    %% Spleen %%
    ODEs(5) = P.Q_Spleen * Y.C_vasc_Artery -  P.Q_Spleen * Y.C_vasc_Spleen;
    %% Heart %%
    ODEs(6) = P.Q_Heart * Y.C_vasc_Artery -  P.Q_Heart * Y.C_vasc_Heart;
    %% Muscle %%
    ODEs(7) = P.Q_Muscle * Y.C_vasc_Artery -  P.Q_Muscle * Y.C_vasc_Muscle;
    %% Adipose %%
    ODEs(8) = P.Q_Adipose * Y.C_vasc_Artery -  P.Q_Adipose * Y.C_vasc_Adipose;


    %% Kidney %%
    renalEliminationPassive = P.GFR * Y.C_plasma_Kidney * P.fu;
    ODEs(9) = P.Q_Kidney * Y.C_vasc_Artery -  P.Q_Kidney * Y.C_vasc_Kidney...
            - renalEliminationPassive;
    
    %% Skin %%
    ODEs(10) = P.Q_Skin * Y.C_vasc_Artery -  P.Q_Skin * Y.C_vasc_Skin;
    %% Brain %%
    ODEs(11) = P.Q_Brain * Y.C_vasc_Artery -  P.Q_Brain * Y.C_vasc_Brain;
    %% Bone %%
    ODEs(12) = P.Q_Bone * Y.C_vasc_Artery -  P.Q_Bone * Y.C_vasc_Bone;

    %% Liver %%
    % input from hepatic artery and portal vein
    liver_input = P.Q_Liver_hep * Y.C_vasc_Artery + P.Q_Gut * Y.C_vasc_Gut ...
                + P.Q_Spleen * Y.C_vasc_Spleen;
    % metabolic clearance
    hepaticMetabolism = P.CLint * Y.C_plasma_Liver * P.fu; %CLu,int in L/s
    
    ODEs(13) = liver_input ...
               - P.Q_Liver * Y.C_vasc_Liver ...
               - hepaticMetabolism;
    
    %% Rest of Body %%
    ODEs(14) = P.Q_RestOfBody * Y.C_vasc_Artery -  P.Q_RestOfBody * Y.C_vasc_RestOfBody;
    
    %% total elimination (for tracking mass conservation) %%
    ODEs(15) = hepaticMetabolism + renalEliminationPassive;
    
    %% convert accumulation rates from moles/s to mg/s
    ODEs(1:15) = ODEs(1:15) * P.RMM * 1e3;  % [mg/s]

end

%%

function Y = calcAlgebraic(m,P,Y)

    % Calculates drug masses and concentrations in all 
    % organ subcompartments at a given time step.

    % m - current total mass of drug in each organ [mg]
    % P - parameter struct
    % Y - struct containing all organ subcompartment concentrations for all
    % organs

    % Cpu - unbound concentration in local plasma [moles/L]
    % C_plasma - total concentration in local plasma [moles/L]
    % m_plasma - total mass in local plasma [mg/L]
    % C_RBC - concentration in local blood cells [moles/L]
    % m_RBC - mass in local blood cells [mg/L]
    % m_vasc - total mass in vascular sub-compartment [m/L]
    % C_vasc - concentration in vascular sub-compartment [moles/L]
    % C_tissue - concentration in tissue subcompartment [moles/L]
    % m_tissue - mass in tissue sub-compartment [moles/L]
    % C_org - concentration in organ [mg/L]

    % Drug mass in each compartment
    mArtery = m(1);
    mVein = m(2);
    mLung = m(3);
    mGut = m(4);
    mSpleen = m(5);
    mHeart = m(6);
    mMuscle = m(7);
    mAdipose = m(8);
    mKidney = m(9);
    mSkin = m(10);
    mBrain = m(11);
    mBone = m(12);
    mLiver = m(13);
    mRestOfBody = m(14);
    
    %% Artery %%   
    Y.Cpu_Artery = Calc_Unbound(mArtery, P.V_vasc_Artery, 0, P.H, P.RMM, ...
                   P.Kpu_RBCs, 0, P.fu); 
    
    Y.C_plasma_Artery = Y.Cpu_Artery / P.fu; 
    
    Y.m_plasma_Artery = Y.C_plasma_Artery * P.V_plasma_Artery * P.RMM * 1e3; 
    Y.C_RBC_Artery = Y.Cpu_Artery * P.Kpu_RBCs; 
    Y.m_RBC_Artery = Y.C_RBC_Artery * P.V_vasc_Artery * P.RMM * 1e3;
    Y.m_vasc_Artery = Y.m_plasma_Artery + Y.m_RBC_Artery; 
    Y.C_vasc_Artery = Y.m_vasc_Artery / P.V_vasc_Artery / P.RMM / 1e3; 
    Y.C_org_Artery = mArtery / P.V_Artery; 
    
    %% Vein %%
    Y.Cpu_Vein = Calc_Unbound(mVein, P.V_vasc_Vein, 0, P.H, P.RMM, ...
                 P.Kpu_RBCs, 0, P.fu);
    
    Y.C_plasma_Vein = Y.Cpu_Vein / P.fu;
    
    Y.m_plasma_Vein = Y.C_plasma_Vein * P.V_plasma_Vein;
    Y.C_RBC_Vein = Y.Cpu_Vein * P.Kpu_RBCs;
    Y.m_RBC_Vein = Y.C_RBC_Vein * P.V_vasc_Vein;
    Y.m_vasc_Vein = Y.m_plasma_Vein + Y.m_RBC_Vein;
    Y.C_vasc_Vein = Y.m_vasc_Vein / P.V_vasc_Vein;
    Y.C_org_Vein = mVein / P.V_Vein;
           
    %% Lung %%
    Y.Cpu_Lung = Calc_Unbound(mLung, P.V_vasc_Lung, P.V_tissue_Lung, P.H, P.RMM, ...
                 P.Kpu_RBCs, P.Kpu_Lung, P.fu);
    
    Y.C_plasma_Lung = Y.Cpu_Lung / P.fu;
    
    Y.m_plasma_Lung = Y.C_plasma_Lung * P.V_plasma_Lung * P.RMM * 1e3;
    Y.C_RBC_Lung = Y.Cpu_Lung * P.Kpu_RBCs;
    Y.m_RBC_Lung = Y.C_RBC_Lung * P.V_vasc_Lung * P.RMM * 1e3;
    Y.m_vasc_Lung = Y.m_plasma_Lung + Y.m_RBC_Lung;
    Y.C_vasc_Lung = Y.m_vasc_Lung / P.V_vasc_Lung / P.RMM / 1e3;  
    Y.C_tissue_Lung = P.Kpu_Lung * Y.Cpu_Lung;
    Y.m_tissue_Lung = Y.C_tissue_Lung * P.V_tissue_Lung * P.RMM * 1e3;
    Y.C_org_Lung = mLung / P.V_Lung;
    
    %% Gut %%
    Y.Cpu_Gut = Calc_Unbound(mGut, P.V_vasc_Gut, P.V_tissue_Gut, P.H, P.RMM, ...
                P.Kpu_RBCs, P.Kpu_Gut, P.fu);
    
    Y.C_plasma_Gut = Y.Cpu_Gut / P.fu;
    
    Y.m_plasma_Gut = Y.C_plasma_Gut * P.V_plasma_Gut * P.RMM * 1e3;
    Y.C_RBC_Gut = Y.Cpu_Gut * P.Kpu_RBCs;
    Y.m_RBC_Gut = Y.C_RBC_Gut * P.V_vasc_Gut * P.RMM * 1e3;
    Y.m_vasc_Gut = Y.m_plasma_Gut + Y.m_RBC_Gut;
    Y.C_vasc_Gut = Y.m_vasc_Gut / P.V_vasc_Gut / P.RMM / 1e3;
    Y.C_tissue_Gut = P.Kpu_Gut * Y.Cpu_Gut;
    Y.m_tissue_Gut = Y.C_tissue_Gut * P.V_tissue_Gut * P.RMM * 1e3;
    Y.C_org_Gut = mGut / P.V_Gut;
    
    %% Spleen %%
    Y.Cpu_Spleen = Calc_Unbound(mSpleen, P.V_vasc_Spleen, ...
                   P.V_tissue_Spleen, P.H, P.RMM, P.Kpu_RBCs, P.Kpu_Spleen, P.fu);
    
    Y.C_plasma_Spleen = Y.Cpu_Spleen / P.fu;
    
    Y.m_plasma_Spleen = Y.C_plasma_Spleen * P.V_plasma_Spleen * P.RMM * 1e3;
    Y.C_RBC_Spleen = Y.Cpu_Spleen * P.Kpu_RBCs;
    Y.m_RBC_Spleen = Y.C_RBC_Spleen * P.V_vasc_Spleen * P.RMM * 1e3;
    Y.m_vasc_Spleen = Y.m_plasma_Spleen + Y.m_RBC_Spleen;
    Y.C_vasc_Spleen = Y.m_vasc_Spleen / P.V_vasc_Spleen / P.RMM / 1e3;  
    Y.C_tissue_Spleen = P.Kpu_Spleen * Y.Cpu_Spleen;
    Y.m_tissue_Spleen = Y.C_tissue_Spleen * P.V_tissue_Spleen * P.RMM * 1e3;
    Y.C_org_Spleen = mSpleen / P.V_Spleen;
    
    %% Heart %%
    Y.Cpu_Heart = Calc_Unbound(mHeart, P.V_vasc_Heart, ...
                    P.V_tissue_Heart, P.H, P.RMM, P.Kpu_RBCs, P.Kpu_Heart, P.fu);
    
    Y.C_plasma_Heart = Y.Cpu_Heart / P.fu;
    
    Y.m_plasma_Heart = Y.C_plasma_Heart * P.V_plasma_Heart * P.RMM * 1e3;
    Y.C_RBC_Heart = Y.Cpu_Heart * P.Kpu_RBCs;
    Y.m_RBC_Heart = Y.C_RBC_Heart * P.V_vasc_Heart * P.RMM * 1e3;
    Y.m_vasc_Heart = Y.m_plasma_Heart + Y.m_RBC_Heart;
    Y.C_vasc_Heart = Y.m_vasc_Heart / P.V_vasc_Heart / P.RMM / 1e3;
    Y.C_tissue_Heart = P.Kpu_Heart * Y.Cpu_Heart;
    Y.m_tissue_Heart = Y.C_tissue_Heart * P.V_tissue_Heart * P.RMM * 1e3;
    Y.C_org_Heart = mHeart / P.V_Heart;
    
    %% Muscle %%
    Y.Cpu_Muscle = Calc_Unbound(mMuscle, P.V_vasc_Muscle, ...
                    P.V_tissue_Muscle, P.H, P.RMM, P.Kpu_RBCs, P.Kpu_Muscle, P.fu);
    
    Y.C_plasma_Muscle = Y.Cpu_Muscle / P.fu;
    
    Y.m_plasma_Muscle = Y.C_plasma_Muscle * P.V_plasma_Muscle * P.RMM * 1e3;
    Y.C_RBC_Muscle = Y.Cpu_Muscle * P.Kpu_RBCs;
    Y.m_RBC_Muscle = Y.C_RBC_Muscle * P.V_vasc_Muscle * P.RMM * 1e3;
    Y.m_vasc_Muscle = Y.m_plasma_Muscle + Y.m_RBC_Muscle;
    Y.C_vasc_Muscle = Y.m_vasc_Muscle / P.V_vasc_Muscle / P.RMM / 1e3;
    Y.C_tissue_Muscle = P.Kpu_Muscle * Y.Cpu_Muscle;
    Y.m_tissue_Muscle = Y.C_tissue_Muscle * P.V_tissue_Muscle * P.RMM * 1e3;
    Y.C_org_Muscle = mMuscle / P.V_Muscle;
    
    %% Adipose %%
    Y.Cpu_Adipose = Calc_Unbound(mAdipose, P.V_vasc_Adipose, ...
                    P.V_tissue_Adipose, P.H, P.RMM, P.Kpu_RBCs, P.Kpu_Adipose, P.fu);
    
    Y.C_plasma_Adipose = Y.Cpu_Adipose / P.fu;
    
    Y.m_plasma_Adipose = Y.C_plasma_Adipose * P.V_plasma_Adipose * P.RMM * 1e3;
    Y.C_RBC_Adipose = Y.Cpu_Adipose * P.Kpu_RBCs;
    Y.m_RBC_Adipose = Y.C_RBC_Adipose * P.V_vasc_Adipose * P.RMM * 1e3;
    Y.m_vasc_Adipose = Y.m_plasma_Adipose + Y.m_RBC_Adipose;
    Y.C_vasc_Adipose = Y.m_vasc_Adipose / P.V_vasc_Adipose / P.RMM / 1e3; 
    Y.C_tissue_Adipose = P.Kpu_Adipose * Y.Cpu_Adipose;
    Y.m_tissue_Adipose = Y.C_tissue_Adipose * P.V_tissue_Adipose * P.RMM * 1e3;
    Y.C_org_Adipose = mAdipose / P.V_Adipose;
    
    %% Kidney %%
    Y.Cpu_Kidney = Calc_Unbound(mKidney, P.V_vasc_Kidney, ...
                    P.V_tissue_Kidney, P.H, P.RMM, P.Kpu_RBCs, P.Kpu_Kidney, P.fu);
    
    Y.C_plasma_Kidney = Y.Cpu_Kidney / P.fu;
    
    Y.m_plasma_Kidney = Y.C_plasma_Kidney * P.V_plasma_Kidney * P.RMM * 1e3;
    Y.C_RBC_Kidney = Y.Cpu_Kidney * P.Kpu_RBCs;
    Y.m_RBC_Kidney = Y.C_RBC_Kidney * P.V_vasc_Kidney * P.RMM * 1e3;
    Y.m_vasc_Kidney = Y.m_plasma_Kidney + Y.m_RBC_Kidney;
    Y.C_vasc_Kidney = Y.m_vasc_Kidney / P.V_vasc_Kidney / P.RMM / 1e3;
    Y.C_tissue_Kidney = P.Kpu_Kidney * Y.Cpu_Kidney;
    Y.m_tissue_Kidney = Y.C_tissue_Kidney * P.V_tissue_Kidney * P.RMM * 1e3;
    Y.C_org_Kidney = mKidney / P.V_Kidney;
    
    %% Skin %%
    Y.Cpu_Skin = Calc_Unbound(mSkin, P.V_vasc_Skin, ...
                    P.V_tissue_Skin, P.H, P.RMM, P.Kpu_RBCs, P.Kpu_Skin, P.fu);
    
    Y.C_plasma_Skin = Y.Cpu_Skin / P.fu;
    
    Y.m_plasma_Skin = Y.C_plasma_Skin * P.V_plasma_Skin * P.RMM * 1e3;
    Y.C_RBC_Skin = Y.Cpu_Skin * P.Kpu_RBCs;
    Y.m_RBC_Skin = Y.C_RBC_Skin * P.V_vasc_Skin * P.RMM * 1e3;
    Y.m_vasc_Skin = Y.m_plasma_Skin + Y.m_RBC_Skin;
    Y.C_vasc_Skin = Y.m_vasc_Skin / P.V_vasc_Skin / P.RMM / 1e3;
    Y.C_tissue_Skin = P.Kpu_Skin * Y.Cpu_Skin;
    Y.m_tissue_Skin = Y.C_tissue_Skin * P.V_tissue_Skin * P.RMM * 1e3;
    Y.C_org_Skin = mSkin / P.V_Skin;
    
    %% Brain %%
    Y.Cpu_Brain = Calc_Unbound(mBrain, P.V_vasc_Brain, ...
                    P.V_tissue_Brain, P.H, P.RMM, P.Kpu_RBCs, P.Kpu_Brain, P.fu);
    
    Y.C_plasma_Brain = Y.Cpu_Brain / P.fu;
    
    Y.m_plasma_Brain = Y.C_plasma_Brain * P.V_plasma_Brain * P.RMM * 1e3;
    Y.C_RBC_Brain = Y.Cpu_Brain * P.Kpu_RBCs;
    Y.m_RBC_Brain = Y.C_RBC_Brain * P.V_vasc_Brain * P.RMM * 1e3;
    Y.m_vasc_Brain = Y.m_plasma_Brain + Y.m_RBC_Brain;
    Y.C_vasc_Brain = Y.m_vasc_Brain / P.V_vasc_Brain / P.RMM / 1e3;
    Y.C_tissue_Brain = P.Kpu_Brain * Y.Cpu_Brain;
    Y.m_tissue_Brain = Y.C_tissue_Brain * P.V_tissue_Brain * P.RMM * 1e3;
    Y.C_org_Brain = mBrain / P.V_Brain;
    
    %% Bone %%
    Y.Cpu_Bone = Calc_Unbound(mBone, P.V_vasc_Bone, ...
                    P.V_tissue_Bone, P.H, P.RMM, P.Kpu_RBCs, P.Kpu_Bone, P.fu);
    
    Y.C_plasma_Bone = Y.Cpu_Bone / P.fu;
    
    Y.m_plasma_Bone = Y.C_plasma_Bone * P.V_plasma_Bone * P.RMM * 1e3;
    Y.C_RBC_Bone = Y.Cpu_Bone * P.Kpu_RBCs;
    Y.m_RBC_Bone = Y.C_RBC_Bone * P.V_vasc_Bone * P.RMM * 1e3;
    Y.m_vasc_Bone = Y.m_plasma_Bone + Y.m_RBC_Bone;
    Y.C_vasc_Bone = Y.m_vasc_Bone / P.V_vasc_Bone / P.RMM / 1e3;
    Y.C_tissue_Bone = P.Kpu_Bone * Y.Cpu_Bone;
    Y.m_tissue_Bone = Y.C_tissue_Bone * P.V_tissue_Bone * P.RMM * 1e3;
    Y.C_org_Bone = mBone / P.V_Bone;
    
    %% Liver %%
    Y.Cpu_Liver = Calc_Unbound(mLiver, P.V_vasc_Liver, ...
                    P.V_tissue_Liver, P.H, P.RMM, P.Kpu_RBCs, P.Kpu_Liver, P.fu);
    
    Y.C_plasma_Liver = Y.Cpu_Liver / P.fu;
    
    Y.m_plasma_Liver = Y.C_plasma_Liver * P.V_plasma_Liver * P.RMM * 1e3;
    Y.C_RBC_Liver = Y.Cpu_Liver * P.Kpu_RBCs;
    Y.m_RBC_Liver = Y.C_RBC_Liver * P.V_vasc_Liver * P.RMM * 1e3;
    Y.m_vasc_Liver = Y.m_plasma_Liver + Y.m_RBC_Liver;
    Y.C_vasc_Liver = Y.m_vasc_Liver / P.V_vasc_Liver / P.RMM / 1e3;
    Y.C_tissue_Liver = P.Kpu_Liver * Y.Cpu_Liver;
    Y.m_tissue_Liver = Y.C_tissue_Liver * P.V_tissue_Liver * P.RMM * 1e3;
    Y.C_org_Liver = mLiver / P.V_Liver; 
    
    %% RestOfBody %%
    Y.Cpu_RestOfBody = Calc_Unbound(mRestOfBody, P.V_vasc_RestOfBody, ...
                    P.V_tissue_RestOfBody, P.H, P.RMM, P.Kpu_RBCs, P.Kpu_RestOfBody, P.fu);
    
    Y.C_plasma_RestOfBody = Y.Cpu_RestOfBody / P.fu;
    
    Y.m_plasma_RestOfBody = Y.C_plasma_RestOfBody * P.V_plasma_RestOfBody * P.RMM * 1e3;
    Y.C_RBC_RestOfBody = Y.Cpu_RestOfBody * P.Kpu_RBCs;
    Y.m_RBC_RestOfBody = Y.C_RBC_RestOfBody * P.V_vasc_RestOfBody * P.RMM * 1e3;
    Y.m_vasc_RestOfBody = Y.m_plasma_RestOfBody + Y.m_RBC_RestOfBody;
    Y.C_vasc_RestOfBody = Y.m_vasc_RestOfBody / P.V_vasc_RestOfBody / P.RMM / 1e3;
    Y.C_tissue_RestOfBody = P.Kpu_RestOfBody * Y.Cpu_RestOfBody;
    Y.m_tissue_RestOfBody = Y.C_tissue_RestOfBody * P.V_tissue_RestOfBody * P.RMM * 1e3;
    Y.C_org_RestOfBody = mRestOfBody / P.V_RestOfBody; 
    
end

%%

function Cu = Calc_Unbound(mass, V_vasc, V_tis, H, RMM, K_RBC, Kpu,fu)

    %Calc_Unbound Calculates unbound concentration in plasma
    %   Employs current drug content of organ (in mg) to determine unbound
    %   concentration (Cu) in blood plasma, assuming equlibrium
    %   distribution between tissue and vascular sub-compartments.
    
    Cu = mass/(V_tis*Kpu + V_vasc*H*K_RBC + V_vasc * (1-H)/fu); %[mg/L]
    Cu = Cu/RMM/1e3; %[moles/L]
    
end
