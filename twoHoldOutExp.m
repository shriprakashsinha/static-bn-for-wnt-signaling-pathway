% function twoHoldOutExp
function twoHoldOutExp(eviDence, model)

% Data Collection
%=====
% Extract data from the gene expression matrix
[uniqueGenes, expressionMatrix, noGenes, noSamples, ...
    groundTruthLabels, transGroundTruthLabels] = readCustomFile('data/geneExpression.mat');

% Load probability values for some of the nodes in the network
fname = dataStorage(model);
load(fname);

% Building the Bayesian Network model
%=====
% Generate directionality between parent and child nodes
[interaction, nodeNames] = generateInteraction(uniqueGenes,model);

% Generate dag for the interaction between nodeNames
dag = mk_adj_mat(interaction, nodeNames, 0);

% To visualise the graphs or bayesian network
gObj = biograph(dag,nodeNames)
gObj = view(gObj);

% BN initialization
N = length(nodeNames); % # of nodes

% Define node sizes. NOTE - nodes are assumed to contain discrete values
nodeSizes = 2*ones(1, N);

% Discrete nodes
discreteNodes = 1:N;

% Create BN
bnet = mk_bnet(dag, nodeSizes, 'names', nodeNames, 'discrete', discreteNodes);

% Hold out expriment
%=====
% Compute no. of positive and negative labels and find indicies of both
noPos = 0;
posLabelIdx = [];
noNeg = 0;
negLabelIdx = [];
for i = 1:noSamples
    if transGroundTruthLabels(i) > 0
        noPos = noPos + 1;
        posLabelIdx = [posLabelIdx, i];
    else
        noNeg = noNeg + 1;
        negLabelIdx = [negLabelIdx, i];
    end
end

% NOTE - Right now using only 2 holdout method
runCnt = 0;
Runs = struct([]);
if ~isempty(strfind(eviDence, 'me'))
    RunsOnObservedMethylation = struct([]);
end

for i = 1:noPos
    for j = 1:noNeg
        % Count for number of runs
        runCnt = runCnt + 1;
        
        % Build test dataset (only 2 examples per test set)
        testDataIdx = [negLabelIdx(j), posLabelIdx(i)];
        dataForTesting = expressionMatrix(:, testDataIdx);
        labelForTesting = transGroundTruthLabels(:, testDataIdx);
        
        % Remove test dataset from the whole dataset and build train
        % dataset
        tmpPosLabelIdx = posLabelIdx;
        tmpNegLabelIdx = negLabelIdx;
        tmpPosLabelIdx(i) = [];
        tmpNegLabelIdx(j) = [];
        trainDataIdx = [tmpNegLabelIdx, tmpPosLabelIdx];
        dataForTraining = expressionMatrix(:, trainDataIdx);
        labelForTraining = transGroundTruthLabels(:, trainDataIdx);
        
        % Define P and CPD for the nodes of the bnet
        cpdStorage = struct([]);
        % Store probabilities for nodes with no parents
        for k = 1:N
            nodeidx = bnet.names(nodeNames{k});
            if isempty(bnet.parents{nodeidx})
                % tables for non-gene measurements
                if ~isempty(strfind(model, 't1'))
                    if strcmp(nodeNames{k}, 'TCF4')
                        cpt = PrTCF4;
                    elseif strcmp(nodeNames{k}, 'MeDKK1')
                        cpt = PrMeDKK1;
                    elseif strcmp(nodeNames{k}, 'MeDACT1')
                        cpt = avgPrMeDACT1;
                    elseif strcmp(nodeNames{k}, 'MeDACT2')
                        cpt = avgPrMeDACT2;
                    elseif strcmp(nodeNames{k}, 'H3k27me3')
                        cpt = avgPrH3K27me3;
                    elseif strcmp(nodeNames{k}, 'H3k4me3')
                        cpt = avgPrH3K4me3;
                    elseif strcmp(nodeNames{k}, 'MeSFRP1')
                        cpt = PrMeSFRP1;
                    elseif strcmp(nodeNames{k}, 'MeSFRP2')
                        cpt = PrMeSFRP2;
                    elseif strcmp(nodeNames{k}, 'MeSFRP4')
                        cpt = PrMeSFRP4;
                    elseif strcmp(nodeNames{k}, 'MeSFRP5')
                        cpt = PrMeSFRP5;
                    elseif strcmp(nodeNames{k}, 'MeWIF1')
                        cpt = PrMeWIF1;
                    elseif strcmp(nodeNames{k}, 'Sample')
                        cpt = PrSample;
                    end
                elseif ~isempty(strfind(model, 't2'))
                    if strcmp(nodeNames{k}, 'TCF4')
                        cpt = PrTCF4;
                    elseif strcmp(nodeNames{k}, 'Sample')
                        cpt = PrSample;
                    end                    
                elseif ~isempty(strfind(model, 'p1'))
                    if strcmp(nodeNames{k}, 'TCF4')
                        cpt = PrTCF4;
                    elseif strcmp(nodeNames{k}, 'BETACAT')
                        cpt = PrBETACAT;
                    end
                end
                cpdStorage(k).parentnode{1} = 0;
                cpdStorage(k).node = nodeNames{k};
                cpdStorage(k).cpt = cpt;
                bnet.CPD{nodeidx} = tabular_CPD(bnet, nodeidx, 'CPT', cpt);
            end
        end
        % Store probabilities for nodes with parents
        gene_cpd = struct([]);
        for k = 1:N
            nodeidx = bnet.names(nodeNames{k});
            if ~isempty(bnet.parents{nodeidx})
                parentidx = bnet.parents{nodeidx};
                noParents = length(parentidx);
                GENE = 0;
                for l = 1:noGenes
                    if strcmp(nodeNames{k}, uniqueGenes{l})
                        % Find cpt of gene parent
                        parent = {};
                        for n = 1:noParents
                            parent{n} = nodeNames{parentidx(n)};
                        end
                        % Assign cpd to parent
                        cnt = 0;
                        parent_cpd = {};
                        for m = 1:length(cpdStorage)
                            for n = 1:noParents
                                if strcmp(parent{n}, cpdStorage(m).node)
                                    cnt = cnt + 1;
                                    parent_cpd{cnt} = cpdStorage(m).cpt;
                                end
                            end
                        end
                        x = generateGenecpd(dataForTraining(l,:),...
                            labelForTraining, uniqueGenes{l}, parent,...
                            parent_cpd, model);
                        gene_cpd(l).vecmedian = x.vecmedian;
                        gene_cpd(l).T = x.T;
                        [r, c] = size(gene_cpd(l).T);
                        cpt = reshape(gene_cpd(l).T, 1, r*c);
                        GENE = 1;
                        break;
                    end
                end
                % tables for non-gene measurements
                if ~GENE
                    if ~isempty(strfind(model, 't1'))
                        if strcmp(nodeNames{k}, 'BETACAT')
                            cpt = PrBETACAT;
                        elseif strcmp(nodeNames{k}, 'DVL2')
                            cpt = PrDVL2;
                        elseif strcmp(nodeNames{k}, 'TRCMPLX')
                            cpt = PrTRCMPLX;
                        end
                    elseif ~isempty(strfind(model, 't2'))
                        if strcmp(nodeNames{k}, 'BETACAT')
                            cpt = PrBETACAT;
                        elseif strcmp(nodeNames{k}, 'DVL2')
                            cpt = PrDVL2;
                        elseif strcmp(nodeNames{k}, 'TRCMPLX')
                            cpt = PrTRCMPLX;
                        end                        
                    elseif ~isempty(strfind(model, 'p1'))
                        if strcmp(nodeNames{k}, 'TRCMPLX')
                            cpt = PrTRCMPLX;
                        end
                    end
                end
                % record the parent index
                for m = 1:noParents
                    cpdStorage(k).parentnode{m} = parentidx(m);
                end
                cpdStorage(k).node = nodeNames{k};
                cpdStorage(k).cpt = cpt;
                bnet.CPD{nodeidx} = tabular_CPD(bnet, nodeidx, 'CPT', cpt);
            end
        end
        % Function to store estimated parameters
        Runs(runCnt).geneCpd = gene_cpd;
        Runs(runCnt).cpdStorage = cpdStorage;
        
        % Function to predict on test data using trained BN
        [r, c] = size(dataForTesting);
        
        switch eviDence
            case 'ge'
                disp(['Testing Example ', num2str(runCnt), ' - Based on all ge']);
                tempTRCMPLXgivenAllge = [];
                if ~isempty(strfind(model, 't'))
                    tempSAMPLE = [];
                end
                
                % Build evidence for inference
                for k = 1:c
                    evidence = cell(1,N);
                    for m = 1:noGenes
                        if dataForTesting(m,k) <= gene_cpd(m).vecmedian
                            evidence{bnet.names(uniqueGenes(m))} = 1;
                        else
                            evidence{bnet.names(uniqueGenes(m))} = 2;
                        end
                    end
                    bnetEngine = jtree_inf_engine(bnet);
                    [engine, loglik] = enter_evidence(bnetEngine, evidence);
                    
                    % Pr(TRCMPLX = 2|gene expression as evidence)
                    margTRCMPLX = marginal_nodes(engine, bnet.names('TRCMPLX'));
                    tempTRCMPLXgivenAllge = [tempTRCMPLXgivenAllge, margTRCMPLX.T(2)];
                    if ~isempty(strfind(model, 't'))
                        % Pr(Sample = 2|gene expression as evidence)
                        margSAMPLE = marginal_nodes(engine, bnet.names('Sample'));
                        tempSAMPLE = [tempSAMPLE, margSAMPLE.T(2)];
                    end
                end
                % Function to store prediction values
                Runs(runCnt).condPrTRCMPLXgivenAllge = tempTRCMPLXgivenAllge;
                if ~isempty(strfind(model,'t'))
                    Runs(runCnt).condPrSAMPLE = tempSAMPLE;
                end
            case 'me'
                disp(['Testing Example ', num2str(runCnt), ' - Based on me']);
                tempTRCMPLX = [];
                tempSAMPLE = [];
                
                % Find the total number of Methylation cases
                methylationIndices = [];
                for k = 1:N
                    if ~isempty(strfind(nodeNames{k}, 'Me'))
                        methylationIndices = [methylationIndices, k];
                    end
                end
                noMethylation = length(methylationIndices);
                
                % Generate binary values for methylation
                q = quantizer( 'mode', 'fixed','roundmode','ceil',...
                    'overflowmode', 'wrap', 'format', [noMethylation 0]);
                noMethylationConfig = 2^noMethylation;
                binaryStatesOfMethylation = {};
                for k = 1:noMethylationConfig
                    binaryStatesOfMethylation{k} = num2bin(q,k);
                end
                
                % Generate methylation evidences
                methylationEvidence = [];
                for k = 1:noMethylationConfig
                    binaryString = binaryStatesOfMethylation{k};
                    for l = 1:noMethylation
                        methylationEvidence(k,l) = str2num(binaryString(l));
                    end
                end
                methylationEvidence = methylationEvidence + 1;
                
                % Build evidence for inference for every test example
                for k = 1:c
                    % For every methylation configuration
                    for l = 1:noMethylationConfig
                        evidence = cell(1,N);
                        % For each of the methylation node
                        cntMe = 0;
                        for m = 1:N
                            if ~isempty(strfind(nodeNames{m}, 'Me'))
                                cntMe = cntMe + 1;
                                evidence{bnet.names(nodeNames(m))} = methylationEvidence(l, cntMe);
                            end
                        end
                        % Build inference engine based on evidence
                        [engine, loglik] = enter_evidence(bnetEngine, evidence);
                        
                        % Pr(TRCMPLX = 2|gene expression as evidence)
                        margTRCMPLX = marginal_nodes(engine, bnet.names('TRCMPLX'));
                        tempTRCMPLX(l,k) =  margTRCMPLX.T(2);
                        if ~isempty(strcmp(model, 't'))
                            % Pr(Sample = 2|gene expression as evidence)
                            margSAMPLE = marginal_nodes(engine, bnet.names('Sample'));
                            tempSAMPLE(l,k) = margSAMPLE.T(2);
                        end
                    end
                end
                % Store predicted results on Observed Methylation
                Runs(runCnt).TRCMPLX = tempTRCMPLX;
                Runs(runCnt).SAMPLE = tempSAMPLE;
        end
    end
end

if strcmp(eviDence, 'ge')
    % Store results
    Results = {};
    cntResult = 0;
    % Estimation of performance levels
    condPrTRCMPLXgivenAllge = [];
    geneEvidence = {};
    if ~isempty(strfind(model, 't'))
        condPrSAMPLE = [];
    end
    labels = [];
    
    % Store prediction values and original labels
    for i = 1:runCnt
        condPrTRCMPLXgivenAllge(i,:) = Runs(i).condPrTRCMPLXgivenAllge;
        geneEvidence{i} = Runs(i).geneEvidence;
        if ~isempty(strfind(model, 't'))
            condPrSAMPLE(i,:) = Runs(i).condPrSAMPLE;
        end
        labels(i,:) = [-1, +1];
    end
    
    % Reshape the vectors
    [r,c] = size(labels);
    labels = reshape(labels, r*c, 1);
    condPrTRCMPLXgivenAllge = reshape(condPrTRCMPLXgivenAllge, r*c, 1);
    if ~isempty(strfind(model, 't'))
        condPrSAMPLE = reshape(condPrSAMPLE, r*c, 1);
    end
    
    % Plot the ROC curve and compute AUC
    [X,Y,T,AUCTRCMPLXgivenAllge] = perfcurve(labels,condPrTRCMPLXgivenAllge,1);
    plot(X,Y,'r');
    xlabel('False positive rate'); ylabel('True positive rate');
    if ~isempty(strfind(model, 't'))
        hold on;
        [X,Y,T,AUCSAMPLE] = perfcurve(labels,condPrSAMPLE,1);
        plot(X,Y,'b');
        legend('TRCMPLX - On', 'SAMPLE - T');
        hold off;
        
        % Perform ks-test the significance between models/evidences/predictions
        [h01, p, ksstat] = kstest2(condPrTRCMPLXgivenAllge, condPrSAMPLE);
    end
    
    if ~isempty(strfind(model, 't1'))
        save('Results.mat', 'Runs', 'condPrTRCMPLXgivenAllge',...
            'condPrNODEgivenIndividualge', 'geneEvidence', 'condPrSAMPLE',...
            'AUCTRCMPLXgivenAllge', 'AUCSAMPLE', 'h01');
    elseif ~isempty(strfind(model, 't2'))
        save('Results.mat', 'Runs', 'condPrTRCMPLXgivenAllge',...
            'condPrNODEgivenIndividualge', 'geneEvidence', 'condPrSAMPLE',...
            'AUCTRCMPLXgivenAllge', 'AUCSAMPLE', 'h01');
    elseif ~isempty(strfind(model, 'p1'))
        save('Results.mat', 'Runs', 'condPrTRCMPLXgivenAllge',...
                'condPrNODEgivenIndividualge', 'geneEvidence', ...
                'AUCTRCMPLXgivenAllge');
    end
else
end