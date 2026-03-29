    % MCResults: PK prediction outputs from MonteCarlo simulations
    % drugParams: biophysicochemical property matrix for molecules
    % preprocessOption: 'pca', 'whiten', or 'raw' for PK data preprocessing
    % featureTransformOption: 'minmax', 'pca','vectorNorm',or 'zscore' 
    %                          transformation of features prior to kmeans clustering
    % clusterFlag: 'y' to perform kmeans clustering. 
    %              'n' to load clusters assignments from file syntheticIDX.mat
    % plotFlag: 'y' to generate analysis plots. 'n' to skip plotting
    
    set(0,'defaulttextinterpreter','latex');
    fontSize = 18; %fontsize for plotting

    % color series to use for figures
    modelColors = {"#cecece","#a559aa","#59a89c","#f0c571","#e02b35","#082a54","#36b700",...
                   "#cecece","#a559aa","#59a89c","#f0c571","#e02b35","#082a54","#36b700"};

    %% names of PBPK models
    modelNames = {'This model', 'Mathew', 'Pearce, calibrated', 'Pearce, uncalibrated'};
    %% model pairing labels for comparitive analyses
    pairNames = {'TM-M', 'TM-PC', 'TM-PU', 'M-PC', 'M-PU', 'PC-PU'};    
    %% parameter names
    ParamLabels = {'RMM', 'pKa,donor', 'pKa,acceptor', 'fu', 'CL$_{int}$',...
                   'logD', 'net charge'};     
    %% Names for PK outcomes of interest
    PKoutcomeNames = {'Vss', 'Cmax,u', 'AUCu', 'T$|$T$>$0'};
    %% ion classes
    ionClassTypes = {'acid', 'base', 'neut', 'zwit'};
    
    %% constants
    Nmodels = length(modelNames);
    Noutcomes = length(PKoutcomeNames);
    Npairings = nchoosek(Nmodels,2); %number of pairwise model comparisons 
    Ndrugs = size(drugParams,1); 

    %% determine molecule charges
    [drugParams, chargeType] = calculateNetCharge(drugParams);

    %% set zeros in the MCResult array to eps
    MCResults(MCResults==0) = eps; 
    
    %% determine number of MonteCarlo realizations
    Nrealizations = size(MCResults,1); 

    %% DATA PRE-PROCESSING
    if strcmp(preprocessOption,'whiten')
        [processedPKpredictions, normalizedPKdata, logTransformedPKdata] = ...
            whiteningPreprocess(MCResults, Noutcomes, Ndrugs, Nrealizations, Nmodels, 'n');
    elseif strcmp(preprocessOption,'pca')
        [processedPKpredictions, normalizedPKdata, logTransformedPKdata] = ...
            PCAPreprocess(MCResults, Noutcomes, Ndrugs, Nrealizations, Nmodels, 'n');
    elseif (strcmp(preprocessOption,'raw'))
        [processedPKpredictions, normalizedPKdata, logTransformedPKdata] = ...
            NormPreprocess(MCResults, Noutcomes, Ndrugs, Nrealizations, Nmodels);
    end


    %% CALCULATE DISSIMILARITY METRICS (Wasserstein distance)
    X = calculateDissimilarityFeatures(processedPKpredictions, Ndrugs, Nmodels, Npairings, Noutcomes);
    
    %% Transform dissimilarity features   
    if strcmp(featureTransformOption,'minmax')
        Xclust = minMaxNorm(X,'col'); %min-max normalize features 
    elseif strcmp(featureTransformOption,'pca') 
        [pcaMat, Xclust, ~, ~, explained,pcaMeans] = pca(X); %PCA-transform features
        save('pcaMat','pcaMat');
        save('pcaMeans','pcaMeans');
    elseif strcmp(featureTransformOption,'vectorNorm')     
        Xclust = X./vecnorm(X); %vector normalize features
    elseif strcmp(featureTransformOption,'zscore') 
        Xclust = zscore(X); %zscore standardize features
    end
    plotFeatureDistributions(Xclust,Noutcomes,Nmodels,PKoutcomeNames,modelNames);

    %% CLUSTER based on model-to-model dissimilarity of Monte Carlo distributions 
    % perform 1000 clustering attempts
    Nattempts = 1000; 
    % allow up to 6 clusters to be assigned in any one clustering attempt
    maxk = 6; 
    % define robust groupings as pairs of molecules that co-locate at least 95% of the time
    frequencyCutoff = 0.95; 

    % execute clustering if desired; alternatively, read in previously
    % saved clusterIDs
    if strcmp(clusterFlag,'y')
        [syntheticIDX, Nclusters, kCentroids, kClusterSizes] = ...
            robustlyCluster(Xclust, Nattempts, maxk, frequencyCutoff, fontSize);  

        save('syntheticIDX','syntheticIDX')
        save('kCentroids','kCentroids')
        save('kClusterSizes','kClusterSizes')
    else
        load syntheticIDX.mat     
        load kCentroids
        load kClusterSizes
        Nclusters = length(unique(syntheticIDX(syntheticIDX>0)));
    end
    
    if strcmp(plotFlag,'y')
        figure(2001);silhouette(Xclust,syntheticIDX);
        %% generate cluster heat map
        % clusterHeatMap(Xclust, syntheticIDX, pairNames, kClusterSizes, fontSize)
        clusterHeatMap(X, syntheticIDX, pairNames, kClusterSizes, fontSize)
    
        %% generate scatterplot of clusters
        %specify desired type of plot (whitened features or dimensionality reduction via pca or tsne)
        % pca, whitened, or tsne
        dimensionFlag = 'pca'; 
        clusterScatterPlot(Xclust, syntheticIDX, dimensionFlag, pairNames, Noutcomes, fontSize);    

        %% perform statistical tests on cluster properties
        moleculeStats(syntheticIDX, drugParams(:,1:6), ParamLabels);

        %% calculate agreement between model predictions and plot
        clusterConcordances(normalizedPKdata, Noutcomes, Nmodels, ...
                            Nclusters, syntheticIDX, PKoutcomeNames, fontSize);   
        
        %% calculate apportionment of chargeTypes amongst clusters
        calculateIonClassApportionment(ionClassTypes, Nclusters, syntheticIDX, ...
                                       drugParams, chargeType, ParamLabels, fontSize);
        
        %% plot predictions for a prototypical molecule from each cluster
        plotClusterPrototypes(Xclust, logTransformedPKdata, kCentroids, syntheticIDX, ...
                              PKoutcomeNames, modelColors, modelNames, fontSize);
    end
    
    %% save arrays for each cluster's parameters
    saveClusterProperties(drugParams(:,1:6), syntheticIDX);

end

%%

function f = minMaxNorm(A,normType)

    if strcmp(normType,'all')
        ub = max(A,[],'all');
        lb = min(A,[],'all');
    elseif strcmp(normType,'col')
        ub = max(A);
        lb = min(A);
    end
    f = (A-lb)./(ub-lb);

end

%%

function [idx, Nclusters, kCentroids, kClusterSizes] = robustlyCluster(Xclust, ...
          Nattempts, maxk, frequencyCutoff, fontSize)

    %% this script performs k-means clustering Nattempt times.
    %% Cluster IDs are then assigned based on a threshold frequency 
    %% (frequencyCutoff) with which samples appear in the same cluster.  
    %% NB: Some samples may not be robustly assigned to any cluster.
    %% these will be assigned cluster ID = 0;

    %% specify clustering statistic
    kmeansApproach = 'silhouette'; 

    %% seed for random number generator
    Npts = size(Xclust,1);
    
    %% Initialize arrays for results
    NclustersAllAttempts = zeros(Nattempts,1);
    adj = zeros(Npts,Npts);
    idxAllAttempts = zeros(Npts,Nattempts);
    pairs = cell(Npts,1);
    idx=zeros(Npts,1);
    
    % perform clustering attempts
    tic
    opts = statset('MaxIter',1000);
    kMeansFunc = @(X,K)(kmeans(X,K,'options',opts,'replicates',5));
    parfor i = 1:Nattempts
        a = load('rngSeeds.mat');
        rng(a.rngSeeds(i));
        eva = evalclusters(Xclust,kMeansFunc,kmeansApproach,...
                          'KList',1:maxk);
        NclustersAllAttempts(i) = eva.OptimalK;
        idxAllAttempts(:,i) = eva.OptimalY; 
    end
    toc
    
    tic
    for i = 1:Nattempts
        for j = 1:Npts-1
            for k = j+1:Npts
                adj(j,k) = adj(j,k) + (idxAllAttempts(j,i) == idxAllAttempts(k,i));
            end
        end
    end
    toc
    adj = adj./Nattempts;
    % save('adj','adj');
    
    %% identify groups based on adjacency frequency cutoff
    %% note that array "grouped" is upper triangular, which avoids duplicate pairings
    grouped = (adj>=frequencyCutoff);
    parfor i = 1:Npts-1
        %identify all molecules that co-cluster with molecule i
        [~,k] = find(grouped(i,:) == 1);        
        pairs{i} = [i * ones(length(k),1),k']; 
    end
    
    %% assign cluster IDs to defined groups
    count=0;
    for i=1:Npts
        if idx(i)==0
            s = unique(pairs{i});
            if ~isempty(s)
                count=count+1;
                idx(s)=count;
            end
        end
    end
    % Identify total number of robust clusters
    Nclusters = max(idx);

    %% reassign cluster IDs in order of decreasing average Wasserstein distance
    [idx, kCentroids, kClusterSizes] = reassignIDs(Nclusters, Xclust, idx, fontSize);

end

%%

function [idx, kCentroids, kClusterSizes] = reassignIDs(Nclusters,...
          Xclust,idx,fontSize)

    kWassersteinMean=zeros(Nclusters,1);
    for i = 1:Nclusters
        kWassersteinMean(i) = mean(Xclust(idx == i,:),'all');
    end
    sortedByDisagreement = flipud(sortrows([kWassersteinMean,(1:Nclusters)'],1));
    sortedByDisagreement = sortedByDisagreement(:,2);
    
    idxNew=zeros(size(idx));
    for i = 1:Nclusters
        idxNew(idx == sortedByDisagreement(i)) = i;
    end
    idx = idxNew;

    %determine number of molecules in each cluster and centroids
    kClusterSizes = zeros(Nclusters,1);
    kCentroids = zeros(Nclusters,size(Xclust,2));
    for i = 1:Nclusters
        kClusterSizes(i) = sum(idx==i); 
        kCentroids(i,:) = mean(Xclust(idx==i,:));        
    end

    %% plot histogram of cluster frequency
    figure(988); histogram(idx,(min(idx):(max(idx)+1))-0.5,'normalization','probability');
    set(gca,'fontsize',fontSize,'xtick',unique(idx),'TickLabelInterpreter','latex');
    xlabel('cluster ID'); xlim([0 max(idx)+1])
    ylabel('frequency');

end

%%

function clusterHeatMap(Xclust, idx, pairNames, kClusterSizes, fontSize)

    reorderedX = sortrows([Xclust,idx],size(Xclust,2)+1);
    reorderedX = reorderedX(:,1:end-1);
    figure(1001);
    imagesc(reorderedX); 
    colormap pink; 
    c = colorbar;
    c.Label.String = 'Wasserstein distance [-]';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 1.2 * fontSize;
    c.TickLabelInterpreter = 'latex';
    c.Label.Rotation = 270;

    groupingIntervals = cumsum(kClusterSizes);
    yline(groupingIntervals,'w','linewidth',2)
    xline(6.5:6:24.5,'w','linewidth',4)

    set(gca,'xtick',[],'ytick',groupingIntervals,'TickLabelInterpreter','latex',...
        'fontsize',fontSize)
    set(gca,'xtick',1:24,'xticklabel',pairNames)    
    ylabel('pseudomolecule ID, sorted by cluster',...
           'Interpreter', 'latex', 'FontSize', 1.2 * fontSize);

end

%%

function clusterScatterPlot(Xclust, idx, dimensionFlag, pairNames, Noutcomes, fontSize)

    Ndims = 6;
    pairs = nchoosek(1:Ndims,2);
    Npairs = size(pairs,1);

    %% select dimensions to plot    
    if strcmp(dimensionFlag,'tsne')
        score = tsne(Xclust,'NumDimensions',Ndims);
                            %,'Perplexity',50,...
                            % 'exaggeration',8, ...
                            % 'distance','mahalanobis');
        chosenAxes = 1:Ndims;
        
    elseif strcmp(dimensionFlag,'pca')
        [~, score, ~, ~, explained] = pca(Xclust);
        chosenAxes = (1:Ndims)';
        chosenVars = explained(1:Ndims);
    else
        score = Xclust;
        sortedScore = flipud(sortrows([(1:(size(Xclust,2)))',...
                             var(Xclust)',...
                             mean(Xclust)'],2));
        chosenAxes = sortedScore(1:Ndims,1);
        chosenVars = sortedScore(1:Ndims,2);
        chosenMeans = sortedScore(1:Ndims,3);
    end
    
    %% make cluster scatter plots
    makeScatterPlots(score,idx,chosenAxes,pairs,fontSize)
    
    %% format plots
    if strcmp(dimensionFlag,'tsne')
        for i=1:Npairs
            nexttile(i);
            xlabel(strcat('dimension$\ $',num2str(pairs(i,1))),'Interpreter','latex')
            ylabel(strcat('dimension$\ $',num2str(pairs(i,2))),'Interpreter','latex')
            set(gca,'TickLabelInterpreter','latex','fontsize',0.7*fontSize)
        end
    elseif strcmp(dimensionFlag,'pca')
        for i=1:Npairs
            nexttile(i);
            xlabel(strcat('PC',num2str(chosenAxes(pairs(i,1))), '$\ $(',...
                   num2str(round(chosenVars(pairs(i,1)))),'$\%$)'),'Interpreter','latex')
            ylabel(strcat('PC',num2str(chosenAxes(pairs(i,2))), '$\ $(',...
                   num2str(round(chosenVars(pairs(i,2)))),'$\%$)'),'Interpreter','latex')
            set(gca,'TickLabelInterpreter','latex','fontsize',0.7*fontSize)
        end
    else
        for i=1:Npairs
            nexttile(i);
            xlabel(strcat('CV$_{',num2str(chosenAxes(pairs(i,1))), '}\ $(',...
                   num2str(round(sqrt(chosenVars(pairs(i,1)))...
                ./ chosenMeans(pairs(i,1)),2)),')'),...
                   'Interpreter','latex')
            ylabel(strcat('CV$_{',num2str(chosenAxes(pairs(i,2))), '}\ $(',...
                   num2str(round(sqrt(chosenVars(pairs(i,2)))...
                ./ chosenMeans(pairs(i,2)),2)),')'),...
                   'Interpreter','latex')
            set(gca,'TickLabelInterpreter','latex','fontsize',0.7*fontSize)
        end
    end

    %% generate errorbarplot of cluster features 
    featureErrorBarPlot(Xclust, idx, Noutcomes, pairNames, fontSize)

end

%%

function agreement = ...
    clusterConcordances(dataMCNorm, Noutcomes, Nmodels, ...
                        Nclusters, idx, PKoutcomeNames, fontSize)

    agreement = zeros(Noutcomes,Nclusters,Nmodels,Nmodels);
    for j = 1:Noutcomes
        y = dataMCNorm{j};   
        for i = 1:Nclusters
            for m = 1:Nmodels
                %select PK data for model m and cluster i (data is log
                %transformed and normalized
                d = squeeze(y(:,m,idx==i));       
                for o = 1:Nmodels
                    %evaluate CCC as measure of prediction agreement across all
                    %molecules and all MC realizations
                    dprime = squeeze(y(:,o,idx==i));
                    agreement(j,i,m,o) = concordanceCorrelation(d(:),dprime(:));
                end
            end
        end 
    end

    %% correlation plots for cluster-specific prediction agreement
    figure(1004); 
    t1004 = tiledlayout(Nclusters, Noutcomes, 'TileSpacing','compact');
    for j = 1:Noutcomes
        for i = 1:Nclusters
            nexttile(tilenum(t1004,i,j));
            imagesc(squeeze(agreement(j,i,:,:)));
            colormap pink;  
            clim([0 1])
            axis square;
            set(gca,'xtick',1:4,'xticklabel',{'T','M','PC','PU'},...
                    'ytick',1:4,'yticklabel',{'T','M','PC','PU'},...
                    'YDir','reverse','XDir','reverse','TickLabelInterpreter','latex',...
                    'fontsize',fontSize);
            if i == 1; title(PKoutcomeNames(j),'fontsize',fontSize); end
            if j == 1; ylabel(strcat('cluster$\ $',num2str(i)),'fontsize',1.2*fontSize); end
        end
    end
    b = colorbar; 
    b.Layout.Tile = 'East';
    b.Label.String = 'concordance correlation';
    b.Label.FontSize = 1.5*fontSize;
    b.Label.Interpreter = 'latex';
    b.TickLabelInterpreter = 'latex';
    b.FontSize = fontSize;
    b.Label.Rotation = 270;

end

%%

function [ionClasses, pChemStats] = calculateIonClassApportionment(ionClassTypes, Nclusters, ...
                                      idx, drugParams, chargeType, ParamLabels, fontSize)

    [Ndrugs,Nparams] = size(drugParams);
    ionClassesInDataset = zeros(length(ionClassTypes),1);
    ionClasses = zeros(length(ionClassTypes),Nclusters,2);    
    
    for i = 1:length(ionClassTypes)
        ionClassesInDataset(i)=sum(strcmp(chargeType(1:Ndrugs),ionClassTypes{i}));
        for j=1:Nclusters
            %fraction of each ion class in the dataset that is each cluster
            ionClasses(i,j,1) = sum(strcmp(chargeType(idx==j),ionClassTypes{i}))...
                              ./ionClassesInDataset(i);
            ionClasses(i,j,2) = sum(strcmp(chargeType(idx==j),ionClassTypes{i}));
        end
    end
    %fraction of cluster that is of each ion class
    ionClasses(:,:,2) = ionClasses(:,:,2)./sum(ionClasses(:,:,2),1); 

    %% plot cluster properties as swarm plots
    pChemStats = plotPropertySwarms(ionClassTypes, ionClasses,...
                 chargeType,drugParams,idx,ParamLabels,fontSize);

end

%%

function plotClusterPrototypes(Xclust, dataMC, kCentroids, idx, ...
         PKoutcomeNames, modelColors, modelNames, fontSize)

    Nclusters = max(idx);
    Noutcomes = length(dataMC);
    Nmodels = size(dataMC{1},2);
    prototypes = zeros(Nclusters,1);
    for i = 1:Nclusters
        d2Centroid = sum((Xclust - kCentroids(i,:)).^2,2);
        d2Centroid(idx~=i) = NaN;
        [~,prototypes(i)] = min(d2Centroid);
    end

    %% plot characteristic molecule distributions for each cluster
    figure(1006); t1006 = tiledlayout(Nclusters,Noutcomes,'tilespacing','compact');
    currentMax = ones(Noutcomes,2) * 0;
    currentMin = ones(Noutcomes,2) * 1e20;
    for i = 1:Noutcomes
        for j = 1:Nclusters
            nexttile(tilenum(t1006,j,i));
            for k = 1:Nmodels
                selectedData = dataMC{i}(:,k,prototypes(j));
                selectedData = selectedData(selectedData>log10(eps));
                nbins = max([floor(range(selectedData))*12,10]);
                h = histogram(selectedData,nbins,...
                           'normalization','pdf',...
                           'displaystyle','stairs','linewidth',2,...
                           'edgecolor',modelColors{k});
                hold on;
    
                %determine axis limits
                % if i == Noutcomes
                %     currentMin(i,1) = min([h.BinEdges(h.BinEdges>-3),currentMin(i,1)]);
                %     currentMax(i,1) = max([h.BinEdges(h.BinEdges<3),currentMax(i,1)]);
                %     currentMax(i,2) = max([h.Values,currentMax(i,2)]);
                % else
                    currentMin(i,1) = min([h.BinEdges,currentMin(i,1)]);
                    currentMax(i,1) = max([h.BinEdges,currentMax(i,1)]);
                    currentMax(i,2) = max([h.Values,currentMax(i,2)]);
                % end
            end
        end

        for j = 1:Nclusters
            nexttile(tilenum(t1006,j,i));    
            xlim([currentMin(i,1) currentMax(i,1)])
            ylabel(strcat('$P$(log',PKoutcomeNames(i),')'))
            xlabel(strcat('log(',PKoutcomeNames(i),')'))
            set(gca,'fontsize',1.2*fontSize,'TickLabelInterpreter','latex')
            axis square;
        end    
    end
    lg=legend(modelNames);
    lg.Interpreter ='latex';
    lg.Layout.Tile = 'North';
    lg.Orientation = 'horizontal';
    lg.FontSize = fontSize*1.2;

end

%%

function [chargetype, netCharge, I] = ionization(pKa_a, pKa_b, pH)

    %% classify whether site should be considered ionizable 
    pKAcidLim = 9.4; %acidic site if pKa_a<9.4
    pKBaseLim = 5.4; %basic site if pKa_b>5.4
    
    %% classify ion class of molecule
    if pKa_a <= pKAcidLim
        chargetype = 'acid';
        if pKa_b > pKBaseLim
            chargetype = 'zwit';
        end
    elseif pKa_b > pKBaseLim
        chargetype = 'base';
    else
        chargetype = 'neut';
    end
    
    %% determine ionization probability for each ionizable site
    ionzn = zeros(1,2);
    pKa = [pKa_a, pKa_b];
    ionizationType = {'acid','base'};
    for i=1:length(pKa)
        ionzn(i) = HendersonHasselbalch(ionizationType{i}, pH, pKa(i))/100; % [%]
    end
    
    %% consider all possible combinations of charged sites and calculate
    %% weighted molecular charge contribution (weighted by probability of each
    %% molecular ionization state)
    Pcharged = ionzn;
    Puncharged = 1-ionzn;
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
    fc = fc_combinations.*Pcharged ...
      + ~fc_combinations.*Puncharged;
    fc = sum(prod(fc,2));

    fa = fa_combinations.*Pcharged ...
      + ~fa_combinations.*Puncharged;
    fa = sum(prod(fa,2));

    fn = fn_combinations.*Pcharged ...
      + ~fn_combinations.*Puncharged;
    fn = sum(prod(fn,2));

    fz = fz_combinations.*Pcharged ...
      + ~fz_combinations.*Puncharged;
    fz = sum(prod(fz,2));
    
    %fraction of molecule population bearing a protonated site
    fwithPos = fwithPos_combinations.*Pcharged ...
            + ~fwithPos_combinations.*Puncharged;
    fwithPos = sum(prod(fwithPos,2));
    
    netCharge = fc-fa;
    
    I.fc = fc;
    I.fa = fa;
    I.fn = fn;
    I.fz = fz;
    I.fwithPos = fwithPos;

end

%%

function f = HendersonHasselbalch(AorB,pH,pKa)
    
    % returns percent of molecules expected to be ionized at site with given
    % pKa in conditions of given pH
    % pKa - acid dissociation constant
    
    if strcmp(AorB,'acid')
        f = 10^(pH-pKa);
    else
        f = 10^(pKa-pH);
    end

    f = f/(1+f)*100;
    
end

%%

function [whitenedForWasserstein, dataMCNorm, dataMC] = whiteningPreprocess(...
          MCResults, Noutcomes, Ndrugs, Nrealizations, Nmodels, existingFlag)
    %% unpack Monte Carlo predictions
    % Nrealizations x Nmodels x Ndrugs
    Vss = (squeeze(MCResults(:,1:4,:)));
    Cmaxu = (squeeze(MCResults(:,5:8,:)));
    AUCu = (squeeze(MCResults(:,9:12,:)));
    T = (squeeze(MCResults(:,13:16,:)));
    
    %% stack outcomes data into a single struct
    dataMC = cell(Noutcomes,1);
    dataMCNorm = cell(Noutcomes,1);
    dataMC{1} = Vss;
    dataMC{2} = Cmaxu;
    dataMC{3} = AUCu;
    dataMC{4} = T;

    %% logtransform and vector normalize PK data 
    for m = 1:Noutcomes 
        dataMC{m} = log10(dataMC{m});
        dataMCNorm{m} = dataMC{m}./norm(dataMC{m}(:));
    end
    
    %% configure data as 2D array for whitening
    rank2Data = zeros(Nmodels*Ndrugs*Nrealizations,Noutcomes);
    modelIndex = zeros(Nmodels*Ndrugs*Nrealizations,1);
    for i = 1:Nmodels
        A = squeeze(dataMC{1}(:,i,:));
        B = squeeze(dataMC{2}(:,i,:));
        C = squeeze(dataMC{3}(:,i,:));
        D = squeeze(dataMC{4}(:,i,:));
        rank2Data((1:Ndrugs*Nrealizations)...
                + (i-1)*Ndrugs*Nrealizations,:) = ...
                  [A(:),B(:),C(:),D(:)];
        modelIndex((1:Ndrugs*Nrealizations)...
                + (i-1)*Ndrugs*Nrealizations) = i;
    end

    %% WHITEN PK DATA to remove correlations
    if strcmp(existingFlag,'y')
        % load existing means and transformation matrix
        load whiteningMeans
        load whiteningMatrix
        whitenedRank2Data = zeros(size(rank2Data));
        for i = 1:Nmodels
            data = rank2Data(modelIndex==i,:);
            rank2Data(modelIndex==i,:) = data - mean(data);
            whitenedRank2Data(modelIndex==i,:) = ...
                rank2Data(modelIndex==i,:) * whiteningMatrix;
            % whitenedRank2Data(modelIndex==i,:) =...
            %    whitenedRank2Data(modelIndex==i,:) + whiteningMeans(i,:);
            whitenedRank2Data(modelIndex==i,:) =...
               whitenedRank2Data(modelIndex==i,:) + mean(data) * whiteningMatrix;
        end
    else
        [whitenedRank2Data, whiteningMeans, whiteningMatrix] = ...
                                pooledVarianceWhitening(rank2Data,modelIndex);
        % [whitenedRank2Data, whiteningMeans, ~, whiteningMatrix] = whiten(rank2Data,1e-8);
        clear rank2Data
        % Preserve means and whitening matrix 
        save('whiteningMeans','whiteningMeans')
        save('whiteningMatrix','whiteningMatrix')
    end
    
    %% reconfigure data to 3D array
    whitenedForWasserstein = zeros(Ndrugs*Nrealizations, Noutcomes, Nmodels);
    for i=1:Nmodels
        whitenedForWasserstein(:,:,i) = ...
            whitenedRank2Data((1:Ndrugs*Nrealizations)+(i-1)*Ndrugs*Nrealizations,:);
    end
    clear whitenedRank2Data
end

%%

function [PCAForWasserstein, dataMCNorm, dataMC] = PCAPreprocess(...
          MCResults, Noutcomes, Ndrugs, Nrealizations, Nmodels, existingFlag)
    %% unpack Monte Carlo predictions
    % Nrealizations x Nmodels x Ndrugs
    Vss = (squeeze(MCResults(:,1:4,:)));
    Cmaxu = (squeeze(MCResults(:,5:8,:)));
    AUCu = (squeeze(MCResults(:,9:12,:)));
    T = (squeeze(MCResults(:,13:16,:)));
    
    %% stack outcomes data into a single struct
    dataMC=cell(Noutcomes,1);
    dataMCNorm=cell(Noutcomes,1);
    dataMC{1} = Vss;
    dataMC{2} = Cmaxu;
    dataMC{3} = AUCu;
    dataMC{4} = T;

    %% logtransform and max-min normalize PK data 
    for m = 1:Noutcomes 
        dataMC{m} = log10(dataMC{m});
        dataMCNorm{m} = dataMC{m}./norm(dataMC{m}(:));
    end
    
    %% configure data as 2D array for whitening
    rank2Data = zeros(Nmodels*Ndrugs*Nrealizations,Noutcomes);
    for i = 1:Nmodels
        A = squeeze(dataMCNorm{1}(:,i,:));
        B = squeeze(dataMCNorm{2}(:,i,:));
        C = squeeze(dataMCNorm{3}(:,i,:));
        D = squeeze(dataMCNorm{4}(:,i,:));
        rank2Data((1:Ndrugs*Nrealizations)...
                + (i-1)*Ndrugs*Nrealizations,:) = ...
                  [A(:),B(:),C(:),D(:)];
    end

    %% WHITEN PK DATA to remove correlations
    if strcmp(existingFlag,'y')
        % load existing means and transformation matrix
        load pcaMeans
        load pcaMatrix
        pcaTransformedRank2Data = (rank2Data - pcaMeans) * pcaMatrix;
    else
        %% matrix for future transformations
        [pcaMatrix, pcaTransformedRank2Data,~,~,explained,pcaMeans] = pca(rank2Data);
        clear rank2Data
        % Preserve means and transformation matrix 
        save('pcaMeans','pcaMeans')
        save('pcaMatrix','pcaMatrix')
    end
    
    %% reconfigure data to 3D array
    PCAForWasserstein = zeros(Ndrugs*Nrealizations, Noutcomes, Nmodels);
    for i=1:Nmodels
        PCAForWasserstein(:,:,i) = ...
            pcaTransformedRank2Data((1:Ndrugs*Nrealizations) + (i-1)*Ndrugs*Nrealizations,:);
    end
    clear pcaTransformedRank2Data
end

%%

function [NormForWasserstein, dataMCNorm, dataMC] = NormPreprocess(...
          MCResults, Noutcomes, Ndrugs, Nrealizations, Nmodels)
    %% unpack Monte Carlo predictions
    % Nrealizations x Nmodels x Ndrugs
    Vss = (squeeze(MCResults(:,1:4,:)));
    Cmaxu = (squeeze(MCResults(:,5:8,:)));
    AUCu = (squeeze(MCResults(:,9:12,:)));
    T = (squeeze(MCResults(:,13:16,:)));
    
    %% stack outcomes data into a single struct
    dataMC=cell(Noutcomes,1);
    dataMCNorm=cell(Noutcomes,1);
    dataMC{1} = Vss;
    dataMC{2} = Cmaxu;
    dataMC{3} = AUCu;
    dataMC{4} = T;

    %% logtransform and max-min normalize PK data 
    for m = 1:Noutcomes 
        dataMC{m} = log10(dataMC{m});
        dataMCNorm{m} = dataMC{m}./norm(dataMC{m}(:)) ; 
    end
    
    %% configure data as 2D array for whitening
    rank2Data = zeros(Nmodels*Ndrugs*Nrealizations,Noutcomes);
    for i = 1:Nmodels
        A = squeeze(dataMCNorm{1}(:,i,:));
        B = squeeze(dataMCNorm{2}(:,i,:));
        C = squeeze(dataMCNorm{3}(:,i,:));
        D = squeeze(dataMCNorm{4}(:,i,:));
        rank2Data((1:Ndrugs*Nrealizations)...
                + (i-1)*Ndrugs*Nrealizations,:) = ...
                  [A(:),B(:),C(:),D(:)];
    end
    
    %% reconfigure data to 3D array
    NormForWasserstein = zeros(Ndrugs*Nrealizations, Noutcomes, Nmodels);
    for i = 1:Nmodels
        NormForWasserstein(:,:,i) = ...
            rank2Data((1:Ndrugs*Nrealizations) + (i-1)*Ndrugs*Nrealizations,:);
    end

end

%%

function X = calculateDissimilarityFeatures(processedPKpredictions, Ndrugs, ...
                                            Nmodels, Npairings, Noutcomes)
    X = zeros(Ndrugs, Npairings*Noutcomes);
    for m = 1:Noutcomes
        for i = 1:Ndrugs
                columnIndex = 0; % column index for distance array
                for j = 1:Nmodels-1
                    for k = j+1:Nmodels
                        columnIndex = columnIndex+1;
                        %% evaluate histogram similarity by Wasserstein distance 
                        X(i,columnIndex + Npairings * (m-1))=...
                            ws_distance(processedPKpredictions((1:1000)+1000*(i-1),m,j),...
                                        processedPKpredictions((1:1000)+1000*(i-1),m,k));
                    end
                end
         end
    end

end

%%

function wsd = ws_distance(u_samples, v_samples, p)
    % WS_DISTANCE 1- and 2- Wasserstein distance between two discrete 
    % probability measures 
    %   
    %   wsd = WS_DISTANCE(u_samples, v_samples) returns the 1-Wasserstein 
    %   distance between the discrete probability measures u and v 
    %   corresponding to the sample vectors u_samples and v_samples
    %
    %   wsd = WS_DISTANCE(u_samples, v_samples, p) returns the p-Wasserstein 
    %   distance between the discrete probability measures u and v
    %   corresponding to the sample vectors u_samples and v_samples. 
    %   p must be 1 or 2.
    %
    % from https://github.com/nklb/wasserstein-distance
    if ~exist('p', 'var')
        p = 1;
    end
    u_samples_sorted = sort(u_samples(:));
    v_samples_sorted = sort(v_samples(:));
    if p == 1
        
        all_samples = unique([u_samples_sorted; v_samples_sorted], 'sorted');
        
        u_cdf = find_interval(u_samples_sorted, all_samples(1:end-1)) ...
            / numel(u_samples);
        v_cdf = find_interval(v_samples_sorted, all_samples(1:end-1)) ...
            / numel(v_samples);
        
        wsd = sum(abs(u_cdf - v_cdf) .* diff(all_samples));
        
    elseif p == 2
        
        u_N = numel(u_samples);
        v_N = numel(v_samples);    
        all_prob = unique([(0:u_N) / u_N, (0:v_N) / v_N], 'sorted').';
        
        u_icdf = u_samples_sorted(fix(all_prob(1:end-1) * u_N) + 1);
        v_icdf = v_samples_sorted(fix(all_prob(1:end-1) * v_N) + 1);
        
        wsd = sqrt(sum((u_icdf-v_icdf).^2 .* diff(all_prob)));
        
    else
        
        error('Only p=1 or p=2 allowed.')
        
    end
end

%%

function idx = find_interval(bounds, vals)
    % Given the two sorted arrays bounds and vals, the function 
    % idx = FIND_INTERVAL(bounds, vals) identifies for each vals(i) the index 
    % idx(i) s.t. bounds(idx(i)) <= vals(i) < bounds(idx(i) + 1).
    m = 0;
    bounds = [bounds(:); inf];
    idx = zeros(numel(vals), 1);
    for i = 1:numel(vals)
        while bounds(m+1) <= vals(i)
            m = m + 1;
        end
        idx(i) = m;
    end
end

%%

function [drugParams, chargeType] = calculateNetCharge(drugParams)
    Ndrugs = size(drugParams,1);
    netCharge = zeros(Ndrugs,1);
    chargeType = cell(Ndrugs,1);
    for i = 1:Ndrugs
        %% determine net charge on each molecule at pH 7.4
        [chargeType{i}, netCharge(i), ~] = ionization(drugParams(i,2),drugParams(i,3),7.4);
    end
    %add charge to parameter array
    drugParams = [drugParams,netCharge]; 
end

%%

function saveClusterProperties(drugParams, idx)

    Nclusters = max(idx);
    for i = 1:Nclusters
        arrayName = strcat('pseudoDrugsPKCluster',num2str(i));
        arrayData = drugParams(idx==i,:);
        eval([strcat(arrayName,'= arrayData;')]);
        save(arrayName,arrayName);
    end

end

%%

% function [Xwh, mu, invMat, whMat] = whiten(X,epsilon)
%     %
%     % Author: Colorado Reed colorado-reed@uiowa.edu
%     if ~exist('epsilon','var')
%         epsilon = 0.0001;
%     end
%     mu = mean(X); 
%     X = bsxfun(@minus, X, mu);
%     A = X'*X;
%     [V,D,~] = svd(A);
%     whMat = sqrt(size(X,1)-1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
%     Xwh = X*whMat;  
%     invMat = pinv(whMat);
% end

%%

function CCC = concordanceCorrelation(A,B)
    sxy=cov([A,B]);
    s2x=sxy(1,1);
    s2y=sxy(2,2);
    sxy=sxy(1,2);
    CCC=2*sxy/(s2x+s2y+(mean(A)-mean(B))^2);
end

%%

function Xclust = GaussRank(X)
    % Gauss Rank
    N = size(X,1);
    Xclust = zeros(size(X));
    for j = 1:size(X,2)
        r = tiedrank(X(:,j));% ranks 1..N
        u = (r - 0.5)/N; % (0,1)
        Xclust(:,j) = sqrt(2)*erfinv(2*u - 1); % ~ N(0,1) marginally
    end
end

%%

function plotFeatureDistributions(Xclust,Noutcomes,Nmodels,PKoutcomeNames,modelNames)
    figure(2000);
    t2000 = tiledlayout(Noutcomes,Nmodels,'tilespacing','compact');
    for i = 1:4 %outcomes
        for j = 1:4 %models
            nexttile(tilenum(t2000,i,j));
            s = j + Noutcomes*(i-1);
            histogram(Xclust(:,s),...
                     'Normalization','pdf','displaystyle','stairs',...
                     'linewidth',1);
            if i == 1
                title(modelNames(j),'interpreter','latex');
            elseif i == Noutcomes
                xlabel(PKoutcomeNames(i),'interpreter','latex');
            end
            set(gca,'fontsize',16,'TickLabelInterpreter','latex')
            ylabel('$P$(y)','interpreter','latex')
        end
    end
end

%%

function pChemStats = plotPropertySwarms(ionClassTypes,ionClasses,...
                      chargeType,drugParams,idx,ParamLabels,fs)

    %constants
    [Ndrugs,Nparams] = size(drugParams);
    Nclusters = length(unique(idx(idx>0)));

    figure(1005); t1005 = tiledlayout(length(ionClassTypes),Nparams,'TileSpacing','compact');
    swarmColors = {"#36b700","#cecece","#a559aa","#59a89c","#f0c571","#e02b35","#082a54",...
                   "#36b700","#cecece","#a559aa","#59a89c","#f0c571","#e02b35","#082a54"};
    swarmMarkers = {'v','o','^','s','d','+','.','x'};

    pChemStats = cell(Nclusters,1);
    for k = 1:Nclusters
        pChemStats{k} = zeros(length(ionClasses),Nparams,3);
        for i = 1:length(ionClassTypes)
            for j = 1:Nparams
                nexttile(tilenum(t1005,i,j));
                selected = and(idx==k,strcmp(chargeType(1:Ndrugs),ionClassTypes{i}));
                p1 = swarmchart(ones(sum(selected),1)*2*k,drugParams(selected,j),24);
                p1.MarkerFaceColor = swarmColors{mod(k,8)+1};
                p1.MarkerEdgeColor = 'k';
                p1.Marker = swarmMarkers(mod(k,8)+1);
                
                hold on;
                xlim([1,2*(Nclusters+2)])
                ylim([min(drugParams(1:Ndrugs,j)),max(drugParams(:,j))]);
                set(gca,'xtick',(1:Nclusters)*2,'xticklabel',[],...
                    'TickLabelInterpreter','latex',...
                    'fontsize',fs)
                if j == 1
                    ylabel(ionClassTypes{i});
                end
                pChemStats{k}(i,j,:) = quantile(drugParams(selected,j),[0.025 0.5 0.975]);
            end    
        end
    end
    
    %% add swarm representing Pearce data
    PearceData = importdata('PearceRatKtpDrugs.csv');
    PearceParams = PearceData.data(:,[1,5,7,3,2,13]);
    NPearceMols = size(PearceParams,1);
    netChargePearce = zeros(NPearceMols,1);
    chargeTypePearce = cell(NPearceMols,1);
    for i = 1:size(PearceParams,1)
        [chargeTypePearce{i},netChargePearce(i),~] = ...
            ionization(PearceParams(i,2),PearceParams(i,3),7.4);
    end
    PearceParams=[PearceParams,netChargePearce];

    for i = 1:4
        for j = 1:Nparams
            nexttile(tilenum(t1005,i,j));
            selected = strcmp(chargeTypePearce,ionClassTypes{i});
            p1 = swarmchart(ones(sum(selected),1)*2*(k+1),PearceParams(selected,j),24);
            p1.MarkerFaceColor = swarmColors{mod(k+1,8)+1};
            p1.MarkerEdgeColor = 'k';
        end
    end

    for j = 1:Nparams
        nexttile(tilenum(t1005,1,j));
        title(ParamLabels(j),'Interpreter','latex');
        nexttile(tilenum(t1005,4,j));
        set(gca,'xtick',(1:Nclusters)*2,'TickLabelInterpreter','latex',...
            'fontsize',fs)
    end
    lg = legend([cellstr(strcat('cluster$\ $',num2str((1:Nclusters)')));cellstr('Pearce dataset')]);
    lg.Interpreter = 'latex';
    lg.Layout.Tile = 'North';
    lg.Orientation = 'horizontal';
    lg.FontSize = 1.2*fs;
    t1005.XLabel.Interpreter = 'latex';
    t1005.XLabel.String = 'cluster';
    t1005.XLabel.FontSize = 1.5*fs;
end

%% 

function featureErrorBarPlot(Xclust, idx, Noutcomes, pairNames, fs)

    Nclusters = length(unique(idx(idx>0)));
    clusterFeatures = zeros(Nclusters,size(Xclust,2),3);
    rows = floor(sqrt(Nclusters));
    cols = ceil(Nclusters/rows);
    figure(1003); 
    t1003 = tiledlayout(rows,cols,'TileSpacing','tight');
    for i = 1:Nclusters
        nexttile(i); 
        clusterFeatures(i,:,:) = quantile(Xclust(idx==i,:),[0.025,0.5,0.975],1)';    
        errorbar(1:size(Xclust,2),clusterFeatures(i,:,2),...
                 clusterFeatures(i,:,2) - clusterFeatures(i,:,1),...
                 clusterFeatures(i,:,3) - clusterFeatures(i,:,2),'o',...
                'linewidth',1)
        for j = 1:Noutcomes
            xline((j-1)*length(pairNames)+0.5,':b');
        end
        set(gca,'fontsize',0.7*fs,'xtick',1:size(Xclust,2),...
           'xticklabel',pairNames,'TickLabelInterpreter','latex');
        ylabel('Wasserstein distance [-]')
    end

end

function makeScatterPlots(score,idx,chosenAxes,pairs,fs)

    figure(1002); 
    Npairs = size(pairs,1);
    r = 3; c = ceil(Npairs/r);
    t1002 = tiledlayout(r,c,"TileSpacing","tight");
    for i = 1:Npairs
        nexttile(i); 
        g = gscatter(score(idx > 0,chosenAxes(pairs(i,1))),...
                     score(idx > 0,chosenAxes(pairs(i,2))),...
                     idx(idx > 0),...
                     [],'....',[8,8,15,8,8]); axis square;
        h = gca; 
        h.Children = circshift(h.Children, -1); %bring cluster 3 to the top of plot
        % g(3).Color='k';
        if i ~= Npairs
            legend('off'); 
        end
    end
    lgd = legend;
    lgd.Layout.Tile = 'North';
    lgd.Interpreter = 'latex';
    lgd.Orientation = 'horizontal';
    lgd.FontSize = fs;

end

%%

function moleculeStats(idx, drugParams, ParamLabels)

    % RGB codes and other formatting for plotting
    markerColors = [[0 0 0];
              [100 100 100];
              [190 190 190];
              [0 0 0];
              [100 100 100];
              [190 190 190]]./256;
    effectSizeLineColor = [67 126 117]./256;
    fs = 18;

    markerChoice = {'o','^','d','s','+','x'};

    Nparams = size(drugParams,2);
    clusters = unique(idx);
    Nclusters = length(clusters);
    Npanels = Nparams * nchoosek(Nclusters,2);
    r = floor(sqrt(Npanels));
    c = ceil(Npanels/r);
    h = zeros(Nclusters,Nclusters,Nparams);
    p = zeros(Nclusters,Nclusters,Nparams);
    ks2stat = zeros(Nclusters,Nclusters,Nparams);
    effect = zeros(Nclusters,Nclusters,Nparams);
    alpha = 0.05/Nparams; % Bonferroni correction for multiple comparisons

    figure(900);
    t900 = tiledlayout(r,c,'tilespacing','tight');
    for k = 1:Nparams
        for i = 1:Nclusters
            X = drugParams(idx == clusters(i),k);
            n1 = length(X);
            for j = i+1:Nclusters
                Y = drugParams(idx == clusters(j),k);
                n2 = length(Y);

                [~, p(i,j,k), ks2stat(i,j,k)] = kstest2(X,Y,'alpha',alpha);
                D = calculateKScriticalValue(alpha,n1,n2);

                h(i,j,k) = (ks2stat(i,j,k) > D); %hypothesis test
                E = meanEffectSize(X,Y,Effect = 'cohen',VarianceType = 'unequal');
                effect(i,j,k) = E.Effect;

                nexttile; 
                L = gardnerAltmanPlot(X,Y,Effect = 'cohen',...
                                      VarianceType = 'unequal');
                axis square;
                % format plot
                for o = 1:2
                    L(o).MarkerEdgeColor = markerColors(o,:);
                    L(o).Marker = markerChoice(o);
                    L(o).LineWidth = 1;
                end
                L(3).Color = effectSizeLineColor;
                L(3).LineWidth = 2;
                
                yyaxis left;
                ylabel(ParamLabels(k),'Interpreter','latex')
                set(gca,'XTickLabel',{num2str(clusters(i)),...
                                      num2str(clusters(j)),"Cohen's d"},...
                                     'TickLabelInterpreter','latex',...
                                     'fontsize',12);
                yyaxis right; set(gca,'YColor',effectSizeLineColor)
                title(''); %remove default title applied by MATLAB
            end
        end
    end
    
    plotEffectSizeBarPlot(clusters,h,effect,ParamLabels,901,'n');
end

%%

function D = calculateKScriticalValue(alpha,n1,n2)
    cAlpha = sqrt(-log(alpha/2)*0.5); % prefactor for KS test critical value
    D = cAlpha * sqrt((n1+n2)/(n1*n2)); % critical value
end

%%

function plotEffectSizeBarPlot(clusters,h,E,ParamLabels,figNum,subgroupFlag)

    d = size(h);
    Nparams = d(end);
    pairings = nchoosek(clusters,2);
    Npairings = size(pairings,1);
    colors = [[235 254 255];
              [8 42 84];
              [68 109 139];
              [144 180 194]]./256;

    rankCriterion = size(E);
    rankCriterion = length(rankCriterion);
    if rankCriterion == 3
        EReshaped = zeros(Npairings,Nparams);
        hReshaped = zeros(Npairings,Nparams);
        for i = 1:Npairings
            for j = 1:Nparams
                EReshaped(i,j) = E(pairings(i,1),pairings(i,2),j);
                hReshaped(i,j) = h(pairings(i,1),pairings(i,2),j);
            end
        end
    else 
        EReshaped = E;
        hReshaped = h;
    end
    
    %% make effect size bar plot
    figure(figNum); hold off;
    b = bar(EReshaped');
    barNumber = size(b,2);
    if  barNumber > 1
        for j = 1:barNumber
            b(j).FaceColor = colors(mod(j,5)+1,:);
        end
    else 
        b.FaceColor = colors(2,:);
    end

    %% Positions for significance indicators
    significanceFlags = hReshaped.*sign(EReshaped).*(abs(EReshaped)+0.2);
    significanceFlags(significanceFlags==0) = NaN;
    hold on;
    pause(5);

    %% Plot significance indicators 
    
    if barNumber > 1
        for i = 1:barNumber 
            x = b(i).XData + b(i).XOffset;
            plot(x,significanceFlags(i,:),'k*','MarkerSize',12);           
        end
    else
        plot(b.XData,...
            significanceFlags(1,:),'k*','MarkerSize',12);
    end

    %% Format plot
    axis square;
    set(gca,'fontsize',18,'TickLabelInterpreter','latex','xticklabel',ParamLabels)
    ylabel("effect size (Cohen's d)",'interpreter','latex','Rotation',90);
    if strcmp(subgroupFlag,'n')
        legend(strcat('cluster$\ $',num2str(clusters(pairings(:,1))),...
              '$\ $vs.$\ $',num2str(clusters(pairings(:,2)))),...
              'interpreter','latex','location','northwest')
    elseif strcmp(subgroupFlag,'y')
        legend(strcat('cluster$\ $',num2str(clusters)),...
              'interpreter','latex','location','northwest')
    end
    xlim([0.25 Nparams+.75])
end

%%

function [Xwh, mu, invMat, whMat] = whiten(X,idx,epsilon)
    %
    % Modified from author: Colorado Reed colorado-reed@uiowa.edu
    if ~exist('epsilon','var')
        epsilon = 0.0001;
    end
    ids = unique(idx)';
    N = zeros(length(ids),1);
    index = 0;
    features = size(X,2);
    sigma = zeros(features,features);

    for i = ids
        index = index + 1;
        N(index) = sum(idx == i);
        A = X(idx==i,:);
        mu = mean(A); 
        A = bsxfun(@minus, A, mu); 
        S = A'*A;
        sigma = sigma + S * (N(index)-1);
    end

    sigma = sigma ./ (sum(N)-length(ids));
    [V,D,~] = svd(sigma);
    whMat = sqrt(size(X,1)-1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
    Xwh = X*whMat;  
    invMat = pinv(whMat);
end

%%

function [whitenedParams, WhiteningMeans, WhiteningMatrix] = pooledVarianceWhitening(data,idx)
    normParams = data./vecnorm(data);
    muParams = zeros(size(normParams));
    for i = unique(idx)'
        muParams(idx==i,:) = muParams(idx==i,:)...
                           + mean(normParams(idx==i,:));
    end
    [whitenedParams, WhMu, ~, WhiteningMatrix] = whiten(normParams - muParams,idx,1e-8); 
    WhMu = (muParams + WhMu) * WhiteningMatrix;
    whitenedParams = whitenedParams + WhMu;

    WhiteningMeans = zeros(length(unique(idx)),size(data,2));
    for i = unique(idx)'
        a = find(idx==i);
        WhiteningMeans(i,:) = WhMu(a(1),:);
    end
end