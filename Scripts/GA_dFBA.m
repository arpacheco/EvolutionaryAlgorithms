function GA_dFBA(saveFileName,dataFile,modelNamesFile,maxNutrients,minGensToConverge,maxGens,popSize,crossoverProb,mutationProb,maxRandIter,objective,objectiveTarget)

% Genetic algorithm heuristic search of environments that confer specific
% microbial community phenotypes. Takes in a .mat structure of dFBA
% phenotypic data and outputs GA results and performance metrics
%
% INPUTS:
%    saveFileName: The name of the .mat file to which to save results
%    dataFile: The dFBA results file. This script uses the following variables within the file:
%       - deltaBiomass: The change in biomass (gDW) in each environment
%       - relAbus: Normalized organism abundances in each environment
%       - alphaDiv: The alpha diversity (species richness) in each environment
%       - shannon: The Shannon entropy of each environment
%       - allMetsFromModels: A list of all metabolites across all simulations
%       - secMetsLog: The secretion flux of each metabolite per organism per simulation
%       - excTable: A table denoting the source organism, the destination organism, and the index of the metabolite (from allMetsFromModels) being exchanged
%    modelNamesFile: A .mat file with the formatted names of genome-scale metabolic models pertaining to the organisms simulated
%    maxNutrients: The maximum number of nutrients per environment combination;
%    minGensToConverge: Minimum number of generations for convergence criterion 3 to be met (no continued improvement)
%    maxGens: The maximum number of GA generations
%    popSize: The population size per generation
%    crossoverProb: Crossover probability ([0:1])
%    mutationProb: Mutation probability  ([0:1])
%    maxRandIter: Number of random seed sets to test
%    objective: An integer from the set [1:5], where each value represents the following optimization targets:
%       - 1: Relative abundance of an organism
%       - 2: Overall Shannon entropy
%       - 3: Total number of metabolic exchanges
%       - 4: Number of exchanges toward an organism
%       - 5: Secretion of a given metabolite
%    objectiveTarget: Necessary only for objectives 1, 4, and 5.
%       - For objective = 1: The index (in modelNamesFile) corresponding to the organism to maximize
%       - For objective = 4: The index (in modelNamesFile) corresponding to the organism to maximize exchange towards
%       - For objective = 5: A string of a metabolite (in allMetsFromModels) whose secretion flux will be optimized
%
% OUTPUTS:
%    Saves all workspace variables in a .mat file
%
% Example, run a GA to identify the environments that maximize the relative abundance of the first organism in the dataset.
%
% GA_dFBA('../Results/GA_dFBA_20NutrientCommunity_MaximizeBSubtilis.mat','../Results/Exhaustive_Community_Parallel_20CS.mat','../Models/modelNamesFacAnaerFormatted.mat',4,10,100,10,0.9,0.35,50,1,1)
%
% Alan R. Pacheco 04/24/2019, last modified 04/28/2021

%% Load data

load(dataFile)
load(modelNamesFile)

% Check for proper target objective
if (objective == 1) || (objective == 4) || (objective == 5)
    if ~exist('objectiveTarget')
        error('Desired objective requires target')
    end     
end

% Define response variable to optimize
responseLabels = {'Relative abundance','Shannon Entropy','Number of exchanges','Number of exchanges toward','Secretion flux of'};
responseLabel = responseLabels{objective};

%% Run genetic algorithm
[meanF,stdF,maxF] = deal(zeros(maxRandIter,maxGens));
gensToConvergence = ones(maxRandIter,1).*maxGens;
maxResponse = zeros(length(objectiveTarget),1);

if strcmp(responseLabel,'Secretion flux of')
    response = sum(squeeze(secMetsLog(:,find(ismember(allMetsFromModels,objectiveTarget)),:)),1)';
elseif strcmp(responseLabel,'Relative abundance')
    response = relAbus(:,objectiveTarget);
elseif strcmp(responseLabel,'Shannon Entropy')
    response = shannon(:,objectiveTarget);
elseif strcmp(responseLabel,'Number of exchanges')
    response = zeros(size(deltaBiomass,1),1);
    for q = 1:size(deltaBiomass,1)
        excTable = excTableLog(:,:,q);
        excTable = excTable(find(sum(excTable,2)),:); % Eliminate padding rows
        excTable = excTable(find(excTable(:,2)-excTable(:,1)),:); % Eliminate self-cycling
        response(q) = size(excTable,1);
    end
elseif contains(responseLabel,'toward')
    response = zeros(size(deltaBiomass,1),1);
    for q = 1:size(deltaBiomass,1)
        excTable = excTableLog(:,:,q);
        excTable = excTable(find(sum(excTable,2)),:); % Eliminate padding rows
        excTable = excTable(find(excTable(:,2)-excTable(:,1)),:); % Eliminate self-cycling
        response(q) = length(find(excTable(:,2) == objectiveTarget));
    end
else
    error('Not a valid response variable type')
end
maxResponse = max(response);

finalCombos = zeros(maxRandIter,1);
for r = 1:maxRandIter

    %Initialize population with random individuals
    P = randi(length(nutrientList),[popSize,maxNutrients]); % P-matrix has all of the individuals (rows: individuals, cols: met indices)
    for rr = 1:size(P,1)
        P(rr,find(rand(1,4) < 0.5)) = 0; % Leads to nutrient combinations of fewer than maxNutrients
        if sum(P(rr,:)) == 0
            P(rr,1) = randi(length(nutrientList),1);
        end
    end
    F = zeros(popSize,1); % F-matrix has the fitness of each of the individuals

    convergenceCount = 0;
    converged = 0;
    for g = 1:maxGens
        comboList = zeros(size(P,1),1);
        for p = 1:popSize
           selectNutrients = unique(P(p,find(P(p,:))));

           % Locate the selected nutrient combination in the results
           sameNumNutrients = find(sum(nutrientMat,2) == length(selectNutrients));
           matches = sameNumNutrients;
            for n = 1:length(selectNutrients)
                matches = intersect(matches,find(nutrientMat(:,selectNutrients(n))));
            end
            selectCombo = matches;
            comboList(p) = selectCombo;

            F(p) = response(selectCombo);
        end

        % Test for convergence
        if g > 1
            [h,~] = ttest2(F,FPrevious,'tail','right');
            if h == 0
                convergenceCount = convergenceCount + 1;
            else
                convergenceCount = 0;
            end
        end
        if convergenceCount >= minGensToConverge
            if abs((max(F)-mean(F))/mean(F)) < mean(F)*0.1... %The difference between the best and the average fitness is less than a given fraction of the fitness of the average individual
                    && abs((max(F) - max((maxF(r,1:g-1))))/max((maxF(r,1:g-1)))) < 0.01 %The difference between the best individual of the current population and the best individual so far is very small
                converged = 1;
            end
        end
        if converged
            gensToConvergence(r) = g-1;
            break
        else

            % Select best individuals as parents
            numParents = ceil(popSize*0.25); % Select top 25% to reproduce (rounded up)
            [~,ranked] = sort(F,'descend');
            parents = P(ranked(1:numParents),:);

            % Create new generation via crossover
            PNew = zeros(popSize,maxNutrients);

            numCarryThrough = floor(popSize*(1-crossoverProb)); % Number of individuals to carry through
            if numCarryThrough > 0
                PNew(1:numCarryThrough,:) = P(ranked(1:numCarryThrough),:);
            end

            parentsCombined = parents(:); % Linearize all mets in parents matrix
            for q = numCarryThrough+1:popSize
                newInd = parentsCombined(randperm(length(parentsCombined),maxNutrients)); % The rest of the individuals are permutations of the parents

                while length(find(newInd)) > length(unique(newInd(find(newInd)))) % If there are duplicates, permute again 
                    newInd = parentsCombined(randperm(length(parentsCombined),maxNutrients));
                end

                PNew(q,:) = newInd;
            end

            % Create variation through mutation
            for i = 1:popSize

                unusedMets = setdiff(1:length(nutrientList),PNew(i,:));

                for j = maxNutrients
                    if rand(1) < mutationProb
                        PNew(i,j) = unusedMets(ceil(rand(1)*length(unusedMets))); % Select randomly from metabolites not used in row
                    end
                end
            end

            % Add random metabolites to empty combinations
            for rr = 1:size(PNew,1)
                if sum(PNew(rr,:)) == 0
                    PNew(rr,1) = randi(length(nutrientList),1);
                end
            end

            % Update P matrix
            P = PNew;

            % Record results
            meanF(r,g) = mean(F);
            stdF(r,g) = std(F);
            maxF(r,g) = max(F);
        end

        FPrevious = F;

    end
    possibleFinalCombos = unique(comboList(find(F == max(F))));
    if length(possibleFinalCombos) > 1
        finalCombos(r) = possibleFinalCombos(1);
    else
        finalCombos(r) = possibleFinalCombos;
    end
end

%% Plot performance
close all
figure
percentiles = prctile(response,[95,99,100]);
percentileColors = [182, 153, 255;35, 61, 202;0, 0, 20]./255;
for i = 1:length(percentiles)
    plot([1:maxGens],ones(maxGens,1).*percentiles(i),'LineWidth',3,'LineStyle','--','Color',percentileColors(i,:))
    hold on
end
for r = 1:maxRandIter
    plot([1:gensToConvergence(r)],maxF(r,1:gensToConvergence(r)),'LineWidth',2,'Color',[135, 179, 141]./255)
    hold on
end
for r = 1:maxRandIter
    if gensToConvergence(r) < maxGens
        scatter(gensToConvergence(r),maxF(r,gensToConvergence(r)),'r')
    end
    hold on
end
legend({'95th percentile','99th percentile','Maximum','GA prediction'})

% Determine stats on when GA passed 99th percentile
passed99th = zeros(size(maxF,1),1);
for i = 1:size(maxF,1)
    passed = min(find(maxF(i,:) > percentiles(2)));
    if ~isempty(passed)
        passed99th(i) = passed;
    else
        passed99th(i) = NaN;
    end
end

disp(['Mean generation passing 99th percentile: ' num2str(nanmean(passed99th)) ' ± ' num2str(nanstd(passed99th)/sqrt(length(passed99th)))])
set(gca,'FontSize',16)
xlabel('Generation')
title('GA Performance')

if contains(responseLabel,'toward')
    ylabel([responseLabel,' ',modelNamesForm{objectiveTarget}])
elseif contains(responseLabel,'of')
    ylabel([responseLabel,' ',objectiveTarget])
else
    ylabel(responseLabel)
end
         
% plot taxonomic compositions of selected combinations
growingOrganisms = find(sum(relAbus(unique(finalCombos),:),1)); % Organisms that appear at least once
colors = parula(length(growingOrganisms));
figure
b = bar(relAbus(unique(finalCombos),growingOrganisms),'stacked','FaceColor','flat');
for k = 1:length(growingOrganisms)
    b(k).CData = colors(k,:);
end
legend(modelNamesForm(growingOrganisms))
ylim([0,1])
xlabel('Nutrient combination')
ylabel('Relative abundance')
set(gca,'FontSize',16)
title('Relative Abundances')

%% Save workspace to a file
save(saveFileName)