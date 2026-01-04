function [Y,Params] = PrepareParameterStructure(Params)

    %% arguments
    % Params - parameter struct
    % Y - struct containing organ subcompartment concentrations for all
    % organs

    %% Read in physiological parameter values for specified species
    if strcmp(Params.species,'hum')
        Physiology = importdata('PhysiologicalParameters.csv');
        Params.bodyWeight = 70; %kg
    elseif strcmp(Params.species,'rat')
        Physiology = importdata('PhysiologicalParametersRat.csv');
        Params.bodyWeight = 0.25; %kg
    end
    Physiology.textdata = Physiology.textdata(2:end,:);

    %% Model compartments and compartment names
    Params.compartmentNumber = Physiology.data(19,1);
    Params.compartmentNames = Physiology.textdata(1:Params.compartmentNumber,1); 
    % add vascular subcompartments
    Params.compartmentNames = [Params.compartmentNames;'RBCs';'Plasma']; 

    %% Organ compartments and names
    Params.tissueNumber = Physiology.data(20,1); 
    Params.tissueNames = Physiology.textdata(3:Params.tissueNumber+2,1);   
    
    %% organ/tissue properties
    density = Physiology.data(1:Params.compartmentNumber,1); %g/cm3
    compartmentVolumes = Physiology.data(1:Params.compartmentNumber,2); % L/kg of bodyweight
    compartmentFlows = Physiology.data(1:Params.compartmentNumber,3); % mL/min/kg of bodyweight^0.75
    CardiacOutput = Physiology.data(21,1); % mL/min/kg of bodyweight^0.75
    vascFraction = Physiology.data(1:Params.compartmentNumber,4); % [-]
    Params.H = Physiology.data(17,4); %[-]

    %% unit conversions and other calculations (e.g., scaling for bodyweight)
    compartmentVolumes = compartmentVolumes...
                       * Params.bodyWeight; % [L]
    compartmentFlows  = compartmentFlows...
                      * Params.bodyWeight^0.75/1e3/60; % [L/s]
    CardiacOutput = CardiacOutput ...
                  * Params.bodyWeight^0.75/1e3/60; % [L/s]
    % Glomerular filtration rate
    Params.GFR = Physiology.data(18,3) ...
               * Params.bodyWeight^0.75/1e3/60; % [L/s]
    Params.PlasmaVol = Physiology.data(16,2) ...
                     * Params.bodyWeight; % [L]
    BloodVolume = Params.PlasmaVol / (1-Params.H); % [L]
    RBCVol = BloodVolume * Params.H; % [L]
    VenousVolume = 0.6 * BloodVolume; % [L]
    ArteryVolume = 0.4 * BloodVolume; % [L]
    compartmentVolumes(1:2) = [VenousVolume; ArteryVolume];
    
    %% Determine tissue partitioning for molecule of interest
    Params = Params.KpuFunc(Params);

    %% Load physiological parameters into parameter structure
    Params.Volumes=[compartmentVolumes; RBCVol]; %[L]    
    for i=1:Params.compartmentNumber
        Params.(strcat('V_',(Params.compartmentNames{i}))) = Params.Volumes(i);
        Params.(strcat('V_vasc_',(Params.compartmentNames{i}))) = ...
                            Params.Volumes(i).*vascFraction(i);
        Params.(strcat('V_tissue_',(Params.compartmentNames{i}))) = ...
                            Params.Volumes(i).*(1-vascFraction(i));
        Params.(strcat('V_plasma_',(Params.compartmentNames{i}))) = ...
                            Params.Volumes(i).*vascFraction(i)*(1-Params.H);
        Params.(strcat('V_RBC_',(Params.compartmentNames{i}))) = ...
                            Params.Volumes(i).*vascFraction(i)*Params.H;
        Params.(strcat('Q_',(Params.compartmentNames{i}))) = compartmentFlows(i);

        %Initial values for drug mass and concetration in each sub-compartment 
        Y.(strcat('C_org_',(Params.compartmentNames{i}))) = 0;
        Y.(strcat('m_vasc_',(Params.compartmentNames{i}))) = 0;
        Y.(strcat('C_vasc_',(Params.compartmentNames{i}))) = 0;
        Y.(strcat('m_plasma_',(Params.compartmentNames{i}))) = 0;
        Y.(strcat('C_plasma_',(Params.compartmentNames{i}))) = 0;
        Y.(strcat('m_RBC_',(Params.compartmentNames{i}))) = 0;
        Y.(strcat('C_RBC_',(Params.compartmentNames{i}))) = 0;
        Y.(strcat('Cpu_',(Params.compartmentNames{i}))) = 0;
    end

    Params.Q_Vein=CardiacOutput;
    Params.Q_Artery=CardiacOutput;
    Params.Q_Lung=CardiacOutput;

    % calculate hepatic artery flowrate
    Params.Q_Liver_hep=Params.Q_Liver-Params.Q_Spleen-Params.Q_Gut; % [L/s]

    %% Estimate properties for 'RestOfBody" compartment
    % RestOfBody represents all organ volume and blood flow not accounted
    % for in the organs that are explicit in the model.
    AvgVascFraction=mean(vascFraction(3:end));
    %calculate blood flows not otherwise accounted for (discrepancy between
    %organs and total cardiac output
    Params.Q_RestOfBody=CardiacOutput...
                -(sum(compartmentFlows(4:Params.compartmentNumber-4))...
                +Params.Q_Liver_hep);
                
    %calculate body mass not otherwise accounted for (discrepancy between
    %organs and total bodyweight)
    m_RestOfBody = Params.bodyWeight - sum(density(3:end) .* compartmentVolumes(3:end)); %kg
    Params.V_RestOfBody = m_RestOfBody; %treating density as ~1g/cm3 
    Params.Volumes(end-1) = Params.V_RestOfBody;

    % assign Kpu, vacular fraction, and IW2EW for RestOfBody as 
    % volume-weighted average over all other tissues (excluding muscle and adipose)
    WeightingTissues = [1:4,7:Params.tissueNumber-1];
    Params.Kpu_RestOfBody = weightedAverage(Params.Volumes(WeightingTissues),...
                                            Params.Kpu(WeightingTissues));
    Params.IW2EW(end-1) = weightedAverage(Params.Volumes(WeightingTissues),...
                                          Params.IW2EW(WeightingTissues));
    AvgVascFraction = weightedAverage(Params.Volumes(WeightingTissues),...
                                      vascFraction(WeightingTissues));

    Params.Kpu(end-1) = Params.Kpu_RestOfBody;
    Params.V_vasc_RestOfBody = ...
                        Params.V_RestOfBody .* AvgVascFraction;
    Params.V_tissue_RestOfBody = ...
                        Params.V_RestOfBody .* (1-AvgVascFraction);
    Params.V_plasma_RestOfBody = ...
                        Params.V_RestOfBody .* AvgVascFraction * (1-Params.H);
    Params.V_RBC_RestOfBody = ...
                        Params.V_RestOfBody .* AvgVascFraction * Params.H;

end

function f = weightedAverage(V,x)
    f = sum(V .* x) / sum(V);
end
