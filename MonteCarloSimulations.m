function MonteCarloSimulations(referenceProperties, selectedTissue, Nsamples, errorScale)
    
    %% arguments
    % referenceProperties - drug-specific properties for molecules of
    %                       interest
    % selectedTissue - numerical ID of tissue of interest
    % Nsamples - number of Monte Carlo realizations to generate
    % errorScale - scaling factor for mean average errors in parameter distributions

    %% specify drug dose, infusion duration, and simulation time
    dose = 15; %[mg/kg]
    dose_duration = mins2secs(1); %[s]
    Tmax = 48; %[hrs] simulation timespan 
    species = 'hum'; % specify species

    Nmolecules = size(referenceProperties,1); %number of molecules to interrogate
    Noutputs = 4; %number of PK outcomes to interrogate (Vss, Cmaxu, AUCu, T0_1uM)
    Nmodels = 4; %number of Kpu models to be evaluated
    
    %% evaluate Monte Carlo outcomes for each molecule
    MCResults = zeros(Nsamples,Noutputs*Nmodels,Nmolecules);
    nominalResults = zeros(Noutputs*Nmodels,Nmolecules);

    for j=1:Nmolecules

        [MCResults(:,:,j),nominalResults(:,j)] = ...
            sampleAndSimulate(referenceProperties(j,:),...
                              selectedTissue, Nsamples,...
                             [dose, dose_duration, Tmax],...
                              errorScale, species);

    end

    save('MCResults','MCResults')
    save('nominalResults','nominalResults')

end

%% 

function [MCoutcomes,nominalOutcomes] = sampleAndSimulate(...
                                        propertyInputs,...
                                        selectedTissue,...
                                        Nsamples, simulationConfig,...
                                        errorScale, species)
    %% arguments
    % propertyInputs - molecule-specific properties
    % Nsamples - desired number of Monte Carlo realizations

    % number of input parameters
    Nparams = length(propertyInputs); 

    %% simulation configuration
    dose = simulationConfig(1);
    dose_duration = simulationConfig(2);
    Tmax = simulationConfig(3);

    %% Type of distribution to sample for each property
    % order: RMM pKa_a pKa_b fu logCLint logD 
    % Beta or Normal 
    sampleType = {'None',...
                  'Beta',...
                  'Beta',...
                  'Beta',...
                  'Norm',...
                  'Norm'};

    %% Upper and lower beta distribution limits [flag as =1 if normal or lognormal]
    maxmin = [[1,14,14,1,   1, 1];...
              [1, 0, 0,1e-4,1, 1]]; 
    
    %% target mean average errors
    MAE = zeros(Nparams,1);
    MAE(2) = 0.42 * errorScale;  %pKa acidic
    MAE(3) = 0.42 * errorScale;  %pKa basic
    MAE(4) = 0.075 * errorScale; %fu
    MAE(5) = 0.25 * errorScale;  %log CLint    
    MAE(6) = 0.467 * errorScale; %logD    

    %% Sample properties from distributions
    distributedParams = zeros(Nsamples,Nparams);
    distributedParams(:,1) = propertyInputs(1); %RMM is fixed
    for i = 2:Nparams
        distributedParams(:,i) = sampleDistribution(propertyInputs(i),...
                                 MAE(i),sampleType{i},Nsamples,maxmin(:,i));
    end
    
    %% convert logCLint to CLint
    distributedParams(:,5) = 10.^distributedParams(:,5); 
    propertyInputs(5) = 10^propertyInputs(5);

    %%  define complete property arrays
    allProperties = [propertyInputs;distributedParams];
    RMM    = allProperties(:,1);
    pKa_a  = allProperties(:,2);
    pKa_b  = allProperties(:,3);    
    fu     = allProperties(:,4);
    CLint  = allProperties(:,5);
    logD74 = allProperties(:,6); 

    %% intialize result arrays for model predictions
    Nmodels = 4;
    Vdeq = zeros(Nsamples+1,Nmodels);
    T0_1uM = zeros(Nsamples+1,Nmodels); 
    Cmaxu = zeros(Nsamples+1,Nmodels);
    AUCu = zeros(Nsamples+1,Nmodels);

    %% run MC simulations 
    parfor i = 1:Nsamples+1
        P=struct;
        P.species = species;
        P.dose = dose;
        P.dose_duration = dose_duration;
        P.RMM = RMM(i);
        P.pKa_a = pKa_a(i);
        P.pKa_b = pKa_b(i);       
        P.fu = fu(i);
        P.CLint = CLint(i);
        P.logD74 = logD74(i);

        %% Run simulations
        P.KpuFunc = @AkpaFarahatKpuModel;
        P1 = SimulatePBPK(P,Tmax); 
        P.KpuFunc = @MathewKpuModel;
        P2 = SimulatePBPK(P,Tmax); 
        P.KpuFunc = @PearceKpuModel;
        P3 = SimulatePBPK(P,Tmax); 
        P.KpuFunc = @PearceNoCorrectionsKpuModel;
        P4 = SimulatePBPK(P,Tmax);

        Vdeq(i,:) = [P1.Vdeq,P2.Vdeq,P3.Vdeq,P4.Vdeq];
        AUCu(i,:) = secs2mins([P1.AUCu(selectedTissue),...
                               P2.AUCu(selectedTissue),...
                               P3.AUCu(selectedTissue),...
                               P4.AUCu(selectedTissue),]); %uM.min
        Cmaxu(i,:) = [P1.Cmaxu(selectedTissue),...
                      P2.Cmaxu(selectedTissue),...
                      P3.Cmaxu(selectedTissue),...
                      P4.Cmaxu(selectedTissue)]; %uM
        T0_1uM(i,:) = [P1.T0_1uM(selectedTissue),...
                     P2.T0_1uM(selectedTissue),...
                     P3.T0_1uM(selectedTissue),...
                     P4.T0_1uM(selectedTissue)]; %hrs    

    end
    
    %% assemble output for all Monte Carlo realizations, all models, and all outcomes 
    MCoutcomes = [Vdeq,Cmaxu,AUCu,T0_1uM];
    nominalOutcomes = MCoutcomes(1,:);
    MCoutcomes = MCoutcomes(2:end,:);

end

%%

function t = secs2mins(x)
    t = x/60;
end

%%

function t = mins2secs(x)
    t = x*60;
end

%%

function distributedParams = sampleDistribution(theta,MAE,sampleType,Nsamples,maxmin)
    if strcmp(sampleType,'Norm')

        distributedParams = generate_dist(theta,MAE,sampleType,Nsamples,maxmin);    

    elseif strcmp(sampleType,'LogN')

        mu = log((theta^2) / sqrt(MAE^2+theta^2));
        sigma = sqrt(log(MAE^2 / (theta^2) + 1));        
        distributedParams = generate_dist(mu,sigma,sampleType,Nsamples,maxmin);

    elseif strcmp(sampleType,'Beta')

        distributedParams = generate_dist(theta,MAE,sampleType,Nsamples,maxmin);

    end
end

%%

function f=generate_dist(mu,MAE,sampleType,N,maxmin)

    if strcmp(sampleType,'Beta')

        % normalize mean and MAE (max-min normalization)
        d = maxmin(1) - maxmin(2);
        mu = (mu - maxmin(2))/d;
        MAE = MAE/d;

        % find beta distribution parameters that provide desired mean and
        % error
        [a,b] = findBeta(mu,MAE);

        % sample from parameterized beta distribution and scale to desired
        % domain
        f = random(sampleType,a,b,[N,size(mu,2)]) * d + maxmin(2);

    else 
        
        sigma = MAE/sqrt(2/pi);
        f = random(sampleType,mu,sigma,[N,size(mu,2)]);

    end

end

%%

function [a,b] = findBeta(mu,sigma)
    
    % find beta distribution parameters a and b that provide a mean
    % and MAE approximating that desired for MC sampling

    test = 1;

    while test > 1e-3

        a = fminbnd(@(a,mu,sigma)assessMAEdiscrepancy(a,mu,sigma),0,1e6,[],mu,sigma);
        b = a * (1-mu) / mu;
        test = assessMAEdiscrepancy(a,mu,sigma);
        
        % if a & b estimates fail to satisfy desired MAE, nudge mean and try again
        if test > 1e-3
            if round(mu) == 0
                mu = mu + 0.01;
            else
                mu = mu - 0.01;
            end
        end

    end
end

%%

function f = assessMAEdiscrepancy(a,mu,MAE)

    %determine whether value of b estimated from a and distribution mean
    %provides the desired MAE for the distribution

    b = a * (1-mu) / mu;
    x = betarnd(a,b,[10000,1]);
    D = mean(abs(x-mu));
    f = abs(D-MAE);

end
