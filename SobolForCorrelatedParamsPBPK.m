function SobolForCorrelatedParamsPBPK(ReferenceData)

    %% argument
    % ReferenceData - array of molecule properties representating the
    % chemical space of interest.  
    % NReference rows (number of molecules) and 6 columns 
    % Columns should be (1) molar mass, (2) donor pKa, (3) acceptor pKa,
    % (4) fraction unbound, (5) logCLint, (6) logD74

    rng('shuffle')

    %% interface function for model to be interrogated
    Model = @PBPK4Sobol;

    %% number of samples for Monte Carlo estimator  
    Nsamples = 2^16; 

    %% number of bootstrapping samples
    NbootstrapSamples = 100;

    %% parameter names
    ParamLabels = {'RMM','pKa$_{donor}$','pKa$_{acceptor}$',...
                   'fu','CLint','logD74','ctrl'};
    % number of parameters
    Nparams = length(ParamLabels); 

    %% outcomes of interest 
    OutcomeLabels = {'Vss','thalf',...
                     'Kpu$_{lung}$','Kpu$_{gut}$','Kpu$_{spleen}$',...
                     'Kpu$_{heart}$','Kpu$_{muscle}$','Kpu$_{adipose}$',...
                     'Kpu$_{kidney}$','Kpu$_{skin}$','Kpu$_{brain}$',...
                     'Kpu$_{bone}$','Kpu$_{liver}$','Kpu$_{RBCs}$',...
                     'Cmaxu$_{lung}$','Cmaxu$_{gut}$','Cmaxu$_{spleen}$',...
                     'Cmaxu$_{heart}$','Cmaxu$_{muscle}$','Cmaxu$_{adipose}$',...
                     'Cmaxu$_{kidney}$','Cmaxu$_{skin}$','Cmaxu$_{brain}$',...
                     'Cmaxu$_{bone}$','Cmaxu$_{liver}$','Cmaxu$_{RBCs}$',...
                     'AUCu$_{lung}$','AUCu$_{gut}$','AUCu$_{spleen}$',...
                     'AUCu$_{heart}$','AUCu$_{muscle}$','AUCu$_{adipose}$',...
                     'AUCu$_{kidney}$','AUCu$_{skin}$','AUCu$_{brain}$',...
                     'AUCu$_{bone}$','AUCu$_{liver}$','AUCu$_{RBCs}$',...
                     'T0_1uM$_{lung}$','T0_1uM$_{gut}$','T0_1uM$_{spleen}$',...
                     'T0_1uM$_{heart}$','T0_1uM$_{muscle}$','T0_1uM$_{adipose}$',...
                     'T0_1uM$_{kidney}$','T0_1uM$_{skin}$','T0_1uM$_{brain}$',...
                     'T0_1uM$_{bone}$','T0_1uM$_{liver}$','T0_1uM$_{RBCs}$',...
                     'ctrl'};
    % number of outcomes
    Noutcomes = length(OutcomeLabels);

    %% initialize matrices for results
    % total correlated indices
    TCj = zeros(Nparams,Noutcomes); 
    % total uncorrelated indices
    TUj = zeros(Nparams,Noutcomes); 
    % total correlated indices for bootstrapped samples
    TCjSub = zeros(NbootstrapSamples,Nparams,Noutcomes); 
    % total correlated indices for bootstrapped samples
    TUjSub = zeros(NbootstrapSamples,Nparams,Noutcomes); 
    %statistical test p-value
    pval = zeros(Nparams,Noutcomes,2); 
    %acceptance/rejection of null hypothesis
    nullHyp = zeros(Nparams,Noutcomes,2); 
    %array of model evaluations
    fmatrix = zeros(Nsamples,Noutcomes,4,Nparams); 

    %% Evaluate covariance matrix for properties 
    %% (add a control parameter, sampled from uniform distribution)
    NReference = size(ReferenceData,1);
    [covMat, Data] = ParameterCopulaRefArray([ReferenceData, rand(NReference,1)]);
    %evaluate parameter means
    mu = mean(Data);

    %% sample from sobolset
    Nparams = size(Data,2);
    p = sobolset(2*Nparams,'skip',Nsamples);
    u = p(1:Nsamples,1:Nparams);
    uPrime = p(1:Nsamples,Nparams+1:end);

    %% calculate indices for all parameters and outcomes   
    for j=1:Nparams
        tic
        % logical selectors for combining samples
        notChosen = ones(1,Nparams); notChosen(j) = 0; 
        notChosen = logical(notChosen);
        chosen = zeros(1,Nparams); chosen(j) = 1; 
        chosen = logical(chosen);

        %Apply selectors to means
        mu_y = mu(chosen);
        mu_z = mu(notChosen);
        
        %apply selectors to covariance matrix
        Sz = covMat(notChosen,notChosen);
        Szy = covMat(notChosen,chosen);
        Sy = covMat(chosen,chosen);
        Syz = covMat(chosen,notChosen);

        %use selectors to split samples
        w = u(:,notChosen);
        vPrime = uPrime(:,chosen);

        %modify sample sets as needed for estimator
        xTilde = norminv(u);
        xTildePrime = norminv(uPrime);
        x = (chol(covMat)'*xTilde')'+mu;
        xPrime = (chol(covMat)'*xTildePrime')'+mu;
        %
        y = x(:,chosen);
        z = x(:,notChosen);
        yPrime = xPrime(:,chosen);
        zPrime = xPrime(:,notChosen); 
        %
        yTildePrime = norminv(vPrime);
        mu_yc = mu_y + Syz*Sz^-1*(z-mu_z)';
        Syc = Sy - Syz*Sz^-1*Szy;
        yBarPrime = (chol(Syc)'*yTildePrime')' + mu_yc';
        %
        zTilde = norminv(w);
        mu_zc1 = mu_z + Szy' * Sy^-1 .* (yPrime - mu_y);
        Szc1 = Sz - Szy * Sy^-1 * Syz;
        z_hat = (chol(Szc1)' * zTilde')' + mu_zc1;          
        %
        yz = zeros(Nsamples, Nparams);
        yz(:, chosen) = y;
        yz(:, notChosen) = z;
        yz = invertCopula(Data,yz,mu,covMat);
        %
        yPrime_zPrime = zeros(Nsamples, Nparams);
        yPrime_zPrime(:, chosen) = yPrime;
        yPrime_zPrime(:, notChosen) = zPrime;
        yPrime_zPrime = invertCopula(Data,yPrime_zPrime,mu,covMat);
        %
        yBarPrimeZ = zeros(Nsamples, Nparams);
        yBarPrimeZ(:, chosen) = yBarPrime;
        yBarPrimeZ(:, notChosen) = z;
        yBarPrimeZ = invertCopula(Data,yBarPrimeZ,mu,covMat);
        %
        yPrimeZHat = zeros(Nsamples, Nparams);
        yPrimeZHat(:, chosen) = yPrime;
        yPrimeZHat(:, notChosen) = z_hat;
        yPrimeZHat = invertCopula(Data,yPrimeZHat,mu,covMat);
        
        %initialize output matrices for parameter j
        f=zeros(Nsamples,Noutcomes,4);

        %solve for function outputs
        parfor n=1:Nsamples 

            %Model returns the outcomes of interest,
            %evaluated for a given parameter set   
            f(n,:,:) = [Model(yz(n,:));...
                        Model(yPrime_zPrime(n,:));...
                        Model(yBarPrimeZ(n,:));...
                        Model(yPrimeZHat(n,:))]';
            
        end

        %determine sensitivity indices for all outcomes of interest
        [TCj(j,:),TUj(j,:)]=calcIndices(f);

        %collect function outputs for bootstrapping
        fmatrix(:,:,:,j)=f;
        toc
    end    

    %% Bootstrapping
    for k=1:NbootstrapSamples
        selected = randsample(Nsamples,Nsamples,true); 
        for j = 1:Nparams 
            [TCjSub(k,j,:),TUjSub(k,j,:)] = calcIndices(squeeze(fmatrix(selected,:,:,j)));
        end
    end
    %evaluate bootstrapping statistics
    confidenceIntervalTCj = quantile(TCjSub,[0.025,0.975],1);
    confidenceIntervalTUj = quantile(TUjSub,[0.025,0.975],1);

    %% ranksum statistical test for equality of medians
    for j=1:Nparams
        for i=1:Noutcomes
            [pval(j,i,1),nullHyp(j,i,1)] = ranksum(TCjSub(:,j,i),TCjSub(:,j,end));
            [pval(j,i,2),nullHyp(j,i,2)] = ranksum(TUjSub(:,j,i),TUjSub(:,j,end));
        end
    end    
     
    %% save Sobol indices and confidence stats
    save('TCj','TCj')
    save('TUj','TUj')
    save('pval','pval')
    save('nullHyp','nullHyp')
    save('confidenceIntervalTCj','confidenceIntervalTCj')
    save('confidenceIntervalTUj','confidenceIntervalTUj')

    %% plot indices
    figure(1); tiledlayout(2,1,'tilespacing','compact');
    nexttile(1); imagesc(TCj.*nullHyp(:,:,1));
    set(gca,'TickLabelInterpreter','latex','fontsize',14,...
            'Xtick',[],...
            'Ytick',1:Nparams,'YTickLabel',ParamLabels);
    c = colorbar; c.Label.Interpreter = 'latex'; c.Label.String = 'TCj';
    c.TickLabelInterpreter = 'latex'; c.Label.Rotation = 270;

    nexttile(2); imagesc(TUj.*nullHyp(:,:,2));
    set(gca,'TickLabelInterpreter','latex','fontsize',14,...
            'Xtick',1:Noutcomes,'XTickLabel',OutcomeLabels,...
            'Ytick',1:Nparams,'YTickLabel',ParamLabels);
    c = colorbar; c.Label.Interpreter = 'latex'; c.Label.String = 'TUj';
    c.TickLabelInterpreter = 'latex'; c.Label.Rotation = 270;
    colormap pink

end

%%

function [TCj,TUj]=calcIndices(f)

    % Calculate the mean
    f0 = mean(f(:,:,1),1);
    % Calculate the variance
    D = (mean(f(:,:,1) .^ 2,1) - f0.^2)+eps;
    
    Sy_56 = f(:,:,2) .* (f(:,:,4) - f(:,:,1)) ./ D;
    Ty_54 = (f(:,:,1) - f(:,:,3)) .^ 2 ./ (2 * D);

    TCj=mean(Sy_56,1);
    TUj=mean(Ty_54,1);

end

%%

function output=PBPK4Sobol(pars)

    P.species = 'hum'; %species ('hum' or 'rat')
    %% drug dose and IV administration duration
    P.dose = 15; %[mg/kg]
    P.dose_duration = mins2sec(1); %[s]
    P.RMM = pars(1);
    P.pKa_a = pars(2);
    P.pKa_b = pars(3);       
    P.fu = pars(4);
    P.CLint = pars(5);    
    P.logD74 = pars(6);
    P.ctrl = pars(7);

    %specify Kpu model
    P.KpuFunc=@AkpaFarahatKpuModel; 

    % execute PBPK simulation
    P=SimulatePBPK(P,48); %48 hr timespan 

    % output outcomes of interest
    output=[P.Vdeq,P.thalf,P.Kp([1:11,13])',P.Cmaxu',P.AUCu',P.T0_1uM',P.ctrl];    
end

%%

function seconds = mins2sec(minutes)
    seconds = minutes * 60;
end

%%

function [covMat,data]=ParameterCopulaRefArray(data)
    rng('shuffle');
    data = data+eps;   

    %% copula fit and sampling
    [Npts,Nprops] = size(data);
    kdata=zeros(Npts,Nprops);
    nudge = 1e-10;
    lb = min(data) - nudge; %lower bound
    ub = max(data) + nudge; %upper bound
    
    %% fit kernel cdf to data
    for i=1:Nprops
        yy=data(:,i);
        [kdata(:,i),~] = ksdensity(yy,yy,'Support',[lb(i),ub(i)],...
                        'BoundaryCorrection','reflection',...
                        'function','cdf');  
    end
    
    %% fit and sample from copula
    kdata(kdata < nudge) = nudge;
    kdata(kdata > 1-nudge) = 1-nudge;
    rhohat = copulafit('Gaussian',kdata);
    x = copularnd('Gaussian',rhohat,10000); 
    covMat=cov(x);

end

%%

function inverted=invertCopula(Data,samples,mu,CovMat)

    % {'RMM','pKa_a','pKa_b','fu','logCLint','logD74'}$; 
    [~,Nprops]=size(Data);  
    nudge = 1e-10;
    lb = min(Data) - nudge; %lower bound
    ub = max(Data) + nudge; %upper bound
    
    %% invert cumulative distribution
    inverted = zeros(size(samples));
    sigmas = sqrt(diag(CovMat))'.*ones(size(samples));
    mus = mu.*ones(size(samples));
    samples = normcdf(samples,mus,sigmas);

    for i=1:Nprops
        [inverted(:,i),~] = ksdensity(Data(:,i),samples(:,i),'Support',[lb(i),ub(i)],...
                           'BoundaryCorrection','reflection',...
                           'function','icdf');  
    end

    inverted(:,5)=10.^inverted(:,5); %convert logCLint to CLint

end
