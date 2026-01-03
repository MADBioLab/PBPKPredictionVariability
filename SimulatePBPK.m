function DrugP = SimulatePBPK(DrugP,Tmax)

    %% arguments
    % DrugP - drug-specific parameters
    % Tmax - simulation timespan of interest, in hours

    %% time span for ODE solution
    tspan =[0,hrs2sec(Tmax)]; %[seconds]
    
    %% Set up structures with physiological parameters etc.
    % Y - struct to contain organ subcompartment concentrations for all
    % organs
    [Y, DrugP] = PrepareParameterStructure(DrugP);
    bodyWeight = DrugP.bodyWeight;
    
    %% initial conditions
    %state variables are masses of drug in each model compartment [mg]
    m0 = zeros(DrugP.compartmentNumber+1,1); %entries 1-13
    
    %% Calculate Vss [L/kg] as Kpu-weighted volume per unit bodyweight    
    % initialize calculation with the weighted contributions from artery
    % and vein compartments (plasma and blood cells)

    volumeContributions = DrugP.PlasmaVol * (1 + DrugP.Kpu_RBCs * DrugP.fu_adj ...
                        * DrugP.H / (1-DrugP.H));

    % iterate over all organs, summing wieghted volume contributions from
    % tissue and vascular sub-compartments
    for i=1:DrugP.tissueNumber
        V_tissue = DrugP.(strcat('V_tissue_',(DrugP.tissueNames{i})));
        V_plasma = DrugP.(strcat('V_plasma_',(DrugP.tissueNames{i})));

        % contribution from tissue sub-compartment [L]
        tissueContribution = DrugP.Kpu(i) * V_tissue * DrugP.fu_adj;

        %contribution from vascular sub-compartment [L]
        vascularContribution = V_plasma * (1 + ...
                             DrugP.Kpu_RBCs * DrugP.fu_adj * DrugP.H ...
                             / (1-DrugP.H)); 
        volumeContributions = volumeContributions ...
                            + tissueContribution ...
                            + vascularContribution;
    end
    DrugP.Vdeq = volumeContributions  / bodyWeight; %[L/kg]

    %% corresponding apparent, steady-state plasma concentration
    DrugP.Css = DrugP.dose/DrugP.Vdeq/DrugP.RMM*1e3; % [uM]
    
    %% set fu to adjusted value and convert Kpu to Kp
    DrugP.fu = DrugP.fu_adj; %[-]
    DrugP.Kp = DrugP.Kpu * DrugP.fu_adj; %[-]
    
    %% Solve system
    ODEs_a = @(time,m,Params,Y,mass)PBPKodes(time,m,Params,Y,mass);
    options = odeset('NonNegative',1:length(m0),'MaxStep',600);
    
    %% solve system of ODEs (state vars in mg)
    [tSim,StateVars] = ode15s(ODEs_a,tspan,m0,options,DrugP,Y,bodyWeight); %% state vars in mg
    
    %% Non compartmental analysis on arterial plasma conc
    arterialPlasmaConc = StateVars(:,1) ...
                      ./ DrugP.V_Artery ...
                      ./ DrugP.B2P; %Arterial plasma concentration in [mg/L]
    [DrugP.Cmax_plasma, ...
     DrugP.AUC_plasma, ...
     DrugP.MRT_plasma,...
     DrugP.CL,DrugP.Vdss,...
     DrugP.thalf] = NCAcalc(arterialPlasmaConc,tSim,DrugP.dose,bodyWeight);
    
    DrugP.C_of_t = [tSim,arterialPlasmaConc./DrugP.RMM*1e3]; % arterial plasma concn in [uM]
  
    CuInterstitialData = cell(DrugP.tissueNumber,1); % conc unbound in intracellular water
    CuIntracellularData = cell(DrugP.tissueNumber,1); % conc unbound in interstitial water
    T1uM = zeros(DrugP.tissueNumber,1); % time above 1 uM
    AUCu = zeros(DrugP.tissueNumber,1); % area under unbound concentration curve
    Cmaxu = zeros(DrugP.tissueNumber,1); % peak unbound concentration curve

    for i = 1:DrugP.tissueNumber
        % vascular sub-compartment
        V_vasc = DrugP.(strcat('V_vasc_',(DrugP.tissueNames{i}))); 
        % tissue sub-compartment
        V_tissue = DrugP.(strcat('V_tissue_',(DrugP.tissueNames{i}))); 
        
        %interstitial unbound concentration [M]
        Cu = Calc_Unbound(StateVars(:,i), ...
                          V_vasc, ...
                          V_tissue, ...
                          DrugP.H, ...
                          DrugP.RMM, ...
                          DrugP.Kpu_RBCs, ...
                          DrugP.Kpu(i), ...
                          DrugP.fu);

        CuInterstitialData{i} = [tSim, Cu*1e6]; %[uM]
        CuIntracellularData{i} = [tSim, Cu*1e6*DrugP.IW2EW(i)]; %[uM]

        % evaluate time above 1 uM
        [j,~] = find(CuIntracellularData{i}(:,2) > 0.1);
        if isempty(j)
            T1uM(i) = eps;
        else
            T1uM(i) = secs2hrs(CuIntracellularData{i}(max(j),1)...
                             - CuIntracellularData{i}(min(j),1)); % [hrs]
        end

        %evaluate AUCu 0-infinity and Cmaxu
        [Cmaxu(i),AUCu(i),~,~,~,~] = NCAcalc(CuIntracellularData{i}(:,2), ...
                                           tSim,DrugP.dose,bodyWeight);
    end
    DrugP.CuInterstitialData = CuInterstitialData;
    DrugP.CuIntracellularData = CuIntracellularData;
    DrugP.T1uM = T1uM; %[hrs]
    DrugP.AUCu = AUCu; %[uM s]
    DrugP.Cmaxu = Cmaxu; %[uM]
    
end

%%

function Cu = Calc_Unbound(mass, V_vasc, V_tis, H, RMM, K_RBC, Kpu,fu)
    %Calc_Unbound Calculates unbound concentration in plasma
    %   Employs current drug content of organ (in mg) to determine unbound
    %   concentration (Cu) in blood plasma, assuming equlibrium
    %   distribution between tissue and vascular sub-compartments.
    
    Cu = mass/(V_tis * Kpu + V_vasc * H * K_RBC + V_vasc * (1-H)/fu); %[mg/L]
    Cu = Cu/RMM/1e3; %[moles/L]

end

%%

function t=hrs2sec(x)
    t=x*3600;
end

%%

function t=secs2hrs(x)
    t=x/3600;
end

%%

function [Cmax,AUC,MRT,CL,Vss,thalf] = NCAcalc(C,t,dose,bodywt)
    %% omit any terminal points that are zeros
    t=t(C>0);C=C(C>0);

    %% estimate terminal decay rate
    % perform linear fit over last 10% of concentration-timecourse
    integ_offset = floor(0.1*length(t));
    pfit = polyfit(t(end-integ_offset:end),log(C(end-integ_offset:end)),1);
    k = -pfit(1); % [1/time]
    thalf = secs2hrs(log(2)/k);
    [Cmax,maxloc]=max(C); %value and location of maximum concentration in array
    
    %% divide data into arrays prior to Cmax and post Cmax
    t1 = t(1:maxloc); C1=C(1:maxloc);
    t2 = t(maxloc+1:end); C2=C(maxloc+1:end);
    
    %% calculate AUC
    % apply normal trapezium rule to integrate pre-Cmax and log version post Cmax
    % and extrapolate to infinity
    Clast = C(end)*exp(-k*t(end));
    AUC = trapz(t1,C1) + trapzlog(t2,C2) + Clast/k; % [conc.time] 
    
    %% calculate MRT (mean residence time)
    y = t.*C(:,1); % [conc.time]
    y1 = y(1:maxloc); 
    y2 = y(maxloc+1:end);
    % first moment of concentration-time
    AUMC = trapz(t1,y1) + trapzlog(t2,y2) + Clast/k * (t(end)+1/k); % [conc.time^2] 
    MRT = AUMC / AUC; % [s]

    %% evaluate systemic clearance and apparent Vss
    CL = dose * bodywt / AUC; % [volume/time]
    Vss = MRT * CL / bodywt; % [volume/kg]

end

%%

function f = trapzlog(x,y)
    f=0;
    dx=diff(x);
    for i=1:length(y)-1
        f = f + dx(i)*(y(i+1)-y(i))/log(y(i+1)/y(i));
    end
end
