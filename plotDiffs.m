function [ fig_handle ] = plotDiffs( mcmcChain, vars_to_plot, param_names )
%% plotDiffs
%   Plot matrix of posterior differences
% INPUT:
%   params
%       is a m x n matrix of parameters with m = number of data points
%       and n = number of parameters
%   paramNames
%       cell array with n strings specifying parameter names
%

% Based on a code that accompanies the book:
% Kruschke, J. K. (2014). Doing Bayesian Data Analysis:
% A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier
% And on the matlab-bayesian-estimation-master package by Nils Winter,
% retrieved from GitHub: https://github.com/NilsWinter/matlab-bayesian-estimation
%-------------------------------------------------------------------------
names = fieldnames(mcmcChain);
if ~exist('param_names','var') 
    param_names = names;
end
% Plot the parameters pairwise, to see correlations:
fig_handle = figure('color','w','NumberTitle','Off','position', [0,0,700,600]);
nVar = numel(vars_to_plot);

for indVar = 1:nVar
    subplot(nVar,nVar, sub2ind([nVar, nVar], indVar, indVar));
    title([param_names{vars_to_plot(indVar)}]);
    for jindVar = 1:nVar
        if jindVar < indVar
            subplot(nVar, nVar, sub2ind([nVar, nVar], indVar, jindVar));
            mbe_plotPost(mcmcChain.(names{vars_to_plot(indVar)})-...
                mcmcChain.(names{vars_to_plot(jindVar)}), 'compVal' ,0);
            title({['\',param_names{vars_to_plot(indVar)}],...
                [' - \',param_names{vars_to_plot(jindVar)}]});
            box off;
            set(gca,'FontSize',6)
        elseif jindVar == indVar
            subplot(nVar,nVar,sub2ind([nVar,nVar],indVar,jindVar));
            ax = subplot(nVar, nVar, sub2ind([nVar, nVar], indVar, jindVar));
            set(ax,'Visible','off');
            text(4,5,['\',param_names{vars_to_plot(indVar)}],...
                'FontWeight','bold','FontSize',14,'HorizontalAlignment','center');
            rectangle('Position',[0 0 10 10]);
%         else
%             subplot(nVar, nVar, sub2ind([nVar, nVar], indVar, jindVar));
%             [r,~] = corr(X(:,indVar), X(:,jindVar));
%             ax = subplot(nVar, nVar, sub2ind([nVar, nVar], indVar, jindVar));
%             rText = num2str(r,'%.3f');
%             text(2,5,['r = ' rText]);
%             set(ax,'visible','off');
%             rectangle('Position',[0 0 10 10]);
        end
    end


end

