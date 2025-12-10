function SobolForCorrelatedParamsPBPK

    rng('shuffle')

    %% model to be interrogated
    Model=@PBPK4Sobol;

    %% number of samples for Monte Carlo estimator  
    Nsamples = 2^16; 

    %% sample sizes for bootstrapping
    subSamples=100;
    subSize=Nsamples;

    %% parameter names
    ParamLabels = {'RMM','pKa_a1','pKa_b1',...
                   'fu','CL','logD74','ctrl'};
    % number of parameters
    Nparams = length(ParamLabels); 

    %% outcomes of interest 
    OutcomeLabels={'VdPlasma','thalf',...
                   'Ktp_{lung}','Ktp_{gut}','Ktp_{spleen}',...
                   'Ktp_{heart}','Ktp_{muscle}','Ktp_{adipose}',...
                   'Ktp_{kidney}','Ktp_{skin}','Ktp_{brain}',...
                   'Ktp_{bone}','Ktp_{liver}','Ktp_{RBCs}',...
                   'Cmaxu_{lung}','Cmaxu_{gut}','Cmaxu_{spleen}',...
                   'Cmaxu_{heart}','Cmaxu_{muscle}','Cmaxu_{adipose}',...
                   'Cmaxu_{kidney}','Cmaxu_{skin}','Cmaxu_{brain}',...
                   'Cmaxu_{bone}','Cmaxu_{liver}','Cmaxu_{RBCs}',...
                   'AUCu_{lung}','AUCu_{gut}','AUCu_{spleen}',...
                   'AUCu_{heart}','AUCu_{muscle}','AUCu_{adipose}',...
                   'AUCu_{kidney}','AUCu_{skin}','AUCu_{brain}',...
                   'AUCu_{bone}','AUCu_{liver}','AUCu_{RBCs}',...
                   'Teff_{lung}','Teff_{gut}','Teff_{spleen}',...
                   'Teff_{heart}','Teff_{muscle}','Teff_{adipose}',...
                   'Teff_{kidney}','Teff_{skin}','Teff_{brain}',...
                   'Teff_{bone}','Teff_{liver}','Teff_{RBCs}',...
                   'control'};
    % number of outcomes
    Noutcomes=length(OutcomeLabels);

    %% initialize matrices for results
    % total correlated indices
    TCj=zeros(Nparams,Noutcomes); 
    % total uncorrelated indices
    TUj=zeros(Nparams,Noutcomes); 
    % total correlated indices for bootstrapped samples
    TCjSub=zeros(subSamples,Nparams,Noutcomes); 
    % total correlated indices for bootstrapped samples
    TUjSub=zeros(subSamples,Nparams,Noutcomes); 
    %statistical test p-value
    pval=zeros(Nparams,Noutcomes,2); 
    %acceptance/rejection of null hypothesis
    nullHyp=zeros(Nparams,Noutcomes,2); 
    %array of model evaluations
    fmatrix=zeros(Nsamples,Noutcomes,4,Nparams); 

    %% Read in molecule properties (model input parameters)
    load pseudoDrugsPKCluster1;
    %% Evaluate covariance matrix for properties 
    %% (add a control parameter, sampled from uniform distribution) 
    [covMat, Data]=ParameterCopulaRefArray([pseudoDrugsPKCluster1,...
                                            rand(size(pseudoDrugsPKCluster1,1),1)]);
    %evaluate parameter means
    mu=mean(Data);

    %% sample from sobolset
    Nparams=size(Data,2);
    p=sobolset(2*Nparams,'skip',Nsamples);
    u=p(1:Nsamples,1:Nparams);
    uPrime=p(1:Nsamples,Nparams+1:end);

    %% calculate indices for all parameters and outcomes   
    for j=1:Nparams
        tic
        % logical selectors for combining samples
        notChosen=ones(1,Nparams); notChosen(j)=0; 
        notChosen=logical(notChosen);
        chosen=zeros(1,Nparams); chosen(j)=1; 
        chosen=logical(chosen);

        %Apply selectors to means
        mu_y=mu(chosen);
        mu_z=mu(notChosen);
        
        %apply selectors to covariance matrix
        Sz=covMat(notChosen,notChosen);
        Szy=covMat(notChosen,chosen);
        Sy=covMat(chosen,chosen);
        Syz=covMat(chosen,notChosen);

        %use selectors to split samples
        w=u(:,notChosen);
        vPrime=uPrime(:,chosen);

        %modify sample sets as needed for estimator
        xTilde=norminv(u);
        xTildePrime=norminv(uPrime);
        x=(chol(covMat)'*xTilde')'+mu;
        xPrime=(chol(covMat)'*xTildePrime')'+mu;
        %
        y=x(:,chosen);
        z=x(:,notChosen);
        yPrime=xPrime(:,chosen);
        zPrime=xPrime(:,notChosen); 
        %
        yTildePrime=norminv(vPrime);
        mu_yc=mu_y + Syz*Sz^-1*(z-mu_z)';
        Syc=Sy - Syz*Sz^-1*Szy;
        yBarPrime=(chol(Syc)'*yTildePrime')' + mu_yc';
        %
        zTilde = norminv(w);
        mu_zc1 = mu_z + Szy' * Sy^-1 .* (yPrime - mu_y);
        Szc1 = Sz - Szy * Sy^-1 * Syz;
        z_hat = (chol(Szc1)' * zTilde')' + mu_zc1;          
        %
        yz = zeros(Nsamples, Nparams);
        yz(:, chosen) = y;
        yz(:, notChosen) = z;
        yz=invertCopula(Data,yz,mu,covMat);
        %
        yPrime_zPrime = zeros(Nsamples, Nparams);
        yPrime_zPrime(:, chosen) = yPrime;
        yPrime_zPrime(:, notChosen) = zPrime;
        yPrime_zPrime=invertCopula(Data,yPrime_zPrime,mu,covMat);
        %
        yBarPrimeZ = zeros(Nsamples, Nparams);
        yBarPrimeZ(:, chosen) = yBarPrime;
        yBarPrimeZ(:, notChosen) = z;
        yBarPrimeZ=invertCopula(Data,yBarPrimeZ,mu,covMat);
        %
        yPrimeZHat = zeros(Nsamples, Nparams);
        yPrimeZHat(:, chosen) = yPrime;
        yPrimeZHat(:, notChosen) = z_hat;
        yPrimeZHat=invertCopula(Data,yPrimeZHat,mu,covMat);
        
        %initialize output matrices for parameter j
        f=zeros(Nsamples,Noutcomes,4);

        %solve for function outputs
        parfor n=1:Nsamples 
            %Model returns the outcomes of interest,
            %evaluated for a given parameter set   
            f(n,:,:)=[Model(yz(n,:));...
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
    for k=1:subSamples
        selected=randsample(Nsamples,subSize,true); 
        for j=1:Nparams 
            [TCjSub(k,j,:),TUjSub(k,j,:)]=calcIndices(squeeze(fmatrix(selected,:,:,j)));
        end
    end
    %evaluate bootstrapping statistics
    confidenceIntervalTCj=quantile(TCjSub,[0.025,0.975],1);
    confidenceIntervalTUj=quantile(TUjSub,[0.025,0.975],1);

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
end

%%

function [Sj,Tj]=calcIndices(f)

    % Calculate the mean
    f0 = mean(f(:,:,1),1);
    % Calculate the variance
    D = (mean(f(:,:,1) .^ 2,1) - f0.^2)+eps;
    
    Sy_56 = f(:,:,2) .* (f(:,:,4) - f(:,:,1)) ./ D;
    Ty_54 = (f(:,:,1) - f(:,:,3)) .^ 2 ./ (2 * D);

    Sj=mean(Sy_56,1);
    Tj=mean(Ty_54,1);

end

%%

function f=PBPK4Sobol(pars)

    P.species='hum'; %species ('hum' or 'rat')
    %% drug dose and IV administration duration
    P.dose = 15; %[mg/kg]
    P.dose_duration = mins2sec(1); %[s]
    P.RMM=pars(1);
    P.pKa_a=pars(2:3);
    P.pKa_b=pars(4:5);       
    P.fu=pars(6);
    P.CLint=pars(7);    
    P.logD74=pars(8);
    P.KpuFunc=@AkpaFarahatKpuModel; %select Kpu model
    P=PBPKSimulation(P); 
    f=[P.Vdeq,P.thalf,P.Kp([1:11,13])',P.Cmaxu',P.AUCu',P.Teff',P.RMM];    
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

    % {'RMM','pKa[a1]','pKa[b1]','fu','CLint','logD74'}; 
    [~,Nprops]=size(Data);  
    nudge=1e-10;
    lb = min(Data) - nudge; %lower bound
    ub = max(Data) + nudge; %upper bound
    
    %% invert cumulative distribution
    inverted=zeros(size(samples));
    sigmas=sqrt(diag(CovMat))'.*ones(size(samples));
    mus=mu.*ones(size(samples));
    samples=normcdf(samples,mus,sigmas);

    for i=1:Nprops
        [inverted(:,i),~]=ksdensity(Data(:,i),samples(:,i),'Support',[lb(i),ub(i)],...
                         'BoundaryCorrection','reflection',...
                         'function','icdf');  
    end

    inverted(:,5)=10.^inverted(:,5); %convert logCLint to CLint

end
