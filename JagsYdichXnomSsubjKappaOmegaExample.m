% JagsYdichXnomSsubjKappaOmegaexample.m
% Based on a code that accompanies the book:
% Kruschke, J. K. (2014). Doing Bayesian Data Analysis:
% A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier
% And on the matlab-bayesian-estimation-master package by Nils Winter,
% retrieved from GitHub: https://github.com/NilsWinter/matlab-bayesian-estimation

% Optional generic preliminaries:
clc; clear all; close all;

%addpath matlab bayesian estimator master, Copyright Nils Winter, 2016
addpath('./mbe');

%output files will be saved in results/output_dir
output_dir ='Jags-YdichXnomSsubjKappaOmega-example';
output_dir = fullfile(pwd,'results',output_dir);
mkdir(output_dir);

%% 1. LOAD DATA
%prepare the data as a struct variable named dataList with fields 'y' (a
%vector with 0 and 1s), 's' (a vector with subject identifiers), 'Ntotal'
%and 'Nsubj'.
myData = readtable(fullfile('data','z6N8z2N7.csv')); % Read data file;
% The y values are in the column named y
y = cellfun(@str2num, myData.y) ; 
% The subject names are in the column named s
s = myData.s;        
[s_labels,~,s]=unique(s); %factor
Ntotal = length(y);  % Compute the total number of flips.
Nsubj = length(s_labels);
dataList = struct('y',y,'s',s,'Ntotal',Ntotal,'Nsubj',Nsubj);


%% 2. SPECIFY MCMC PROPERTIES

% Number of MCMC steps that are saved for EACH chain
numSavedSteps = 3334;

% Number of separate MCMC chains
nChains = 3;

% Number of steps that are thinned, matjags will only keep every nth
% step. This does not affect the number of saved steps. I.e. in order
% to compute 10000 saved steps, matjags/JAGS will compute 50000 steps
% If memory isn't an issue, Kruschke recommends to use longer chains
% and no thinning at all.
thinSteps = 1;

% Number of burn-in samples
burnInSteps = 1000;

% The parameters that are to be monitored
parameters = {'theta','kappa','omega'};

% Initial values of MCMC chains based on data:
for i=1:length(s_labels)
    theta(i) = mean(y(s==i));
end
omega = mean(y);
kappa = 100;

% Set initial values for latent variable in each chain
for i=1:nChains
    initsList(i) = struct('theta', theta, 'omega',omega,'kappaMinusTwo',kappa-2);
end


%% 3. SPECIFY THE JAGS MODEL
% This will write a JAGS model to a text file
% You can also write the JAGS model directly to a text filemodelString = ...

modelString = {...
    'model {'                                                               ,...
    'for ( i in 1:Ntotal ) {'                                               ,...
    '   y[i] ~ dbern( theta[s[i]] )'                                        ,...
    '  }'                                                                   ,...
    'for (sIdx in 1:Nsubj) {'                                               ,...
    '  theta[sIdx] ~ dbeta( (kappa-2)*omega + 1 , (kappa-2)*(1-omega)+1 )'  ,...
    '  }'                                                                   ,...
    'kappa <- kappaMinusTwo+2'                                              ,...
    'kappaMinusTwo ~ dgamma(0.01,0.01)'                                     ,...
    'omega ~ dbeta(1,1)'                                                    ,...
    '}'                                                                     };

% close quote for modelString
fid = fopen(fullfile(output_dir,'model.txt'),'w'); fprintf(fid, '%s\r\n',modelString{:});
fclose(fid); model = fullfile(output_dir,'model.txt');


%% 4. RUN THE CHAINS USING MATJAGS
% In case you have the Parallel Computing Toolbox, use ('doParallel',1)
[~, ~, mcmcChain] = matjags(...
    dataList,...
    model,...
    initsList,...
    'monitorparams', parameters,...
    'nChains', nChains,...
    'nBurnin', burnInSteps,...
    'thin', thinSteps,...
    'verbosity',1,...
    'nSamples',numSavedSteps);

% This transforms the output of matjags into the format that mbe is
% using
mcmcChain = mbe_restructChains(mcmcChain);
%Change fieldnames to meaningful parameter names. parameters should be
%named in the order they are listed in the 'parameters' variable
mcmcChain = struct2cell(mcmcChain);
mcmcChain = cell2struct(mcmcChain, [strcat('theta',s_labels)','kappa','omega']);
%% 5. EXAMINE THE CHAINS
mbe_diagMCMC(mcmcChain);

%% 6. EXAMINE THE RESULTS
% At this point, we want to use all the chains at once, so we
% need to concatenate the individual chains to one long chain first
mcmcChain = mbe_concChains(mcmcChain);
% Get parameter names
params = fieldnames(mcmcChain);
% Get summary and posterior plots of each parameter;
for i = 1:numel(params)
    summary.(params{i}) = mbe_summary(mcmcChain.(params{i}));
end
%plot data pairs and posterior distributions
fig_handle = mbe_plotPairs(mcmcChain,numSavedSteps);
saveas(fig_handle,fullfile(output_dir,'pairs.png'));

% Plot posterior distribution for each parameter separately
for i = 1:numel(params)
        fig_handle = figure('color','w','NumberTitle','Off','position', [0,0,700,600]);    
        mbe_plotPost(mcmcChain.(params{i}),'credMass',0.95,...
        'xlab',['\',params{i}]);
        saveas(fig_handle,[params{i},'.png']);
end

%plot differences between selected variables
fig_handle = plotDiffs(mcmcChain,[1,2]);
saveas(fig_handle, fullfile(output_dir,'diffs.png'));

%% SAVE FILES
save(fullfile(output_dir,'workspace.mat'),'summary','mcmcChain')
