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
output_dir ='BattingAverage';
output_dir = fullfile(pwd,'results',output_dir);
mkdir(output_dir);

%% 1. LOAD DATA
%prepare the data as a struct variable named dataList with fields 'z' (a
%vector with number of successes), 'N' (a vector with number of flips), 
%'s' (a vector with numeric subject identifiers), 'c' (a vector with
%numeric category identifiers), 'Nc' (number of categories) and 'Ns' (number
%of subjects)
myData = readtable(fullfile('data','BattingAverage.csv')); % Read data file;
z = myData.Hits ; 
N = myData.AtBats;
% The subject names are in the column named s
s_labels = myData.Player;        
% s_labels = cellfun(@(str)str(str~=' '), s_labels, 'UniformOutput',0);
s = reshape(1:(length(s_labels)), size(s_labels));
c = myData.PriPos;
[c_labels,~,c]=unique(c); %factor
% % %remove spaces from primary position labels
% % c_labels = cellfun(@(str)str(str~=' '), c_labels, 'UniformOutput',0);

Ns = length(s_labels);  % Compute the total number of subjects
Nc = length(c_labels);  % Compute the number of categories

dataList = struct('z',z,'N',N,'sub',s,'cat',c,'Nsub',Ns,'Ncat',Nc);


%% 2. SPECIFY MCMC PROPERTIES

% Number of MCMC steps that are saved for EACH chain
numSavedSteps = 334;

% Number of separate MCMC chains
nChains = 3;

% Number of steps that are thinned, matjags will only keep every nth
% step. This does not affect the number of saved steps. I.e. in order
% to compute 10000 saved steps, matjags/JAGS will compute 50000 steps
% If memory isn't an issue, Kruschke recommends to use longer chains
% and no thinning at all.
thinSteps = 10;

% Number of burn-in samples
burnInSteps = 100;

% The parameters that are to be monitored
parameters = {'theta','kappa','omega','kappa0','omega0'};
parameter_names = [strcat('theta ',s_labels)',...
    strcat('kappa ', c_labels)', strcat('omega ',c_labels)', 'kappa0','omega0'];

% Initial values of MCMC chains based on data:
for i=1:Ns
    theta(i) = 0.01+z(i)/N(i)*0.98;
end
for i=1:Nc
    omega(i) = 0.01+(sum(z(c==i))/sum(N(c==i)))*0.98;
    kappa(i) = 100;                             %chosen randomly
end
omega0 = 0.01+sum(z)/sum(N)*0.98;
kappa0 = 100;

% Set initial values for latent variable in each chain
for i=1:nChains
    initsList(i) = struct('theta', theta, 'omega',omega,'kappaMinusTwo',...
        kappa-2,'omega0',omega0,'kappaMinusTwo0',kappa0-2);
end


%% 3. SPECIFY THE JAGS MODEL
% This will write a JAGS model to a text file
% You can also write the JAGS model directly to a text filemodelString = ...

modelString = {...
    'model {'                                                               ,...
    'for ( sIdx in 1:Nsub ) {'                                                   ,...
    '   z[sIdx] ~ dbin( theta[sIdx], N[sIdx] )'                                      ,...
    '   theta[sIdx] ~ dbeta(omega[cat[sIdx]]*(kappa[cat[sIdx]]-2)+1,'                    ,...
    '                    (1-omega[cat[sIdx]])*(kappa[cat[sIdx]]-2)+1)'                 ,...
    '  }'                                                                   ,...
    'for (cIdx in 1:Ncat) {'                                                     ,...
    '   omega[cIdx] ~ dbeta(omega0*(kappa0 - 2)+1,'                              ,...
    '                    (1-omega0)*(kappa0-2)+1)'                           ,...
    '   kappa[cIdx] = kappaMinusTwo[cIdx]+2'                                      ,...
    '   kappaMinusTwo[cIdx] ~ dgamma(0.01,0.01) # mean=1 , sd=10 (generic vague)',...
    '  }'                                                                   ,...
    'omega0 ~ dbeta(1,1)'                                                   ,...
    'kappa0 = kappaMinusTwo0 + 2'                                           ,...
    'kappaMinusTwo0 ~ dgamma(0.01, 0.01) #generic vague'                    ,...
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

%% 5. EXAMINE THE CHAINS
% mbe_diagMCMC(mcmcChain);

%% 6. EXAMINE THE RESULTS
% At this point, we want to use all the chains at once, so we
% need to concatenate the individual chains to one long chain first
mcmcChain = mbe_concChains(mcmcChain);
% Get parameter names
%error if number of monitored parameters and number of parameter names do
%not match
params = fieldnames(mcmcChain);
if numel(parameter_names)~=numel(fieldnames(mcmcChain))
    error(sprintf('number of monitored parameters is %d while number of parameter names is %d',...
        numel(fieldnames(mcmcChain)),numel(parameter_names)));
end
% Get summary and posterior plots of each parameter;
for i = 1:numel(parameter_names)
    summary.(params{i}) = mbe_summary(mcmcChain.(params{i}));
end
%plot data pairs and posterior distributions
fig_handle = mbe_plotPairs(mcmcChain,numSavedSteps,...
    [1,2,3,length(params)-1,length(params)],parameter_names);
saveas(fig_handle,fullfile(output_dir,'pairsPlot.png'));

% Plot posterior distribution for each parameter separately
for i = [1,10,numel(params)-10,numel(params)] %chosen randomly
        fig_handle = figure('color','w','NumberTitle','Off','position', [0,0,700,600]);    
        mbe_plotPost(mcmcChain.(params{i}),'credMass',0.95,...
        'xlab',['\',parameter_names{i}]);
        saveas(fig_handle,fullfile(output_dir,[parameter_names{i},'.png']));
end

%plot differences between selected variables
fig_handle = plotDiffs(mcmcChain,[1,2],parameter_names);
saveas(fig_handle, fullfile(output_dir,'diffs.png'));

%% SAVE FILES
save(fullfile(output_dir,'workspace.mat'),'summary','mcmcChain')
