clc; clear all; close all;
%addpath matlab bayesian estimator master, Copyright Nils Winter, 2016
addpath('./mbe');

fileNameRoot='BernMetrop';                       % for output filenames

%Specify the data, to be used in the likelihood function.
myData = [repelem(0,6),repelem(1,14)];

% Define the relative probability of the target distribution,
% as a function of vector theta. For our application, this
% target distribution is the unnormalized posterior distribution.
targetRelProb = @( theta , data ) (computeLikelihood( theta , data )...
    .* computePrior( theta ));

%Specify the length of the trajectory, i.e., the number of jumps to try:
trajLength = 5000; % arbitrary large number
%Initialize the vector that will store the results:
trajectory = nan(1 , trajLength);
%Specify where to start the trajectory:
trajectory(1) = 0.01; % arbitrary value
%Specify the burn-in period:
burnIn = ceil( 0.0 * trajLength ); % arbitrary number, less than trajLength
% Initialize accepted, rejected counters, just to monitor performance:
nAccepted = 0;
nRejected = 0;

% Now generate the random walk. The 't' index is time or trial in the walk.
% Specify seed to reproduce same random walk:
rng(47405);
% Specify standard deviation of proposal distribution:
proposalSD = [0.02,0.2,2.0];
proposalSD = proposalSD(2);
for t = 1:trajLength-1
    currentPosition = trajectory(t);
    % Use the proposal distribution to generate a proposed jump.
    proposedJump = normrnd( 0 , proposalSD );
    % Compute the probability of accepting the proposed jump.
    probAccept = min( 1, targetRelProb( currentPosition + proposedJump , myData )...
        / targetRelProb( currentPosition , myData ) );
    % Generate a random uniform value from the interval [0,1] to
    % decide whether or not to accept the proposed jump.
    if rand < probAccept
        %accept the proposed jump
        trajectory(t+1) = currentPosition + proposedJump;
        %increment the accepted counter, just to monitor performance
        if t > burnIn; nAccepted = nAccepted + 1; end
    else
        % reject the proposed jump, stay at current position
        trajectory(t+1) = currentPosition;
        % increment the rejected counter, just to monitor performance
        if  t > burnIn;  nRejected = nRejected + 1; end
    end
end

% Extract the post-burnIn portion of the trajectory.
acceptedTraj = trajectory((burnIn+1) : length(trajectory));

% End of Metropolis algorithm.

%-----------------------------------------------------------------------
% Display the chain
figure('NumberTitle','Off','Color','w','Units', 'Centimeters', 'Position', [0,0,8,16]);

% Posterior histogram:
subplot(3,1,1);
[~,neff,~,~,~,~,~] = psrf(acceptedTraj');
paramInfo = mbe_plotPost( acceptedTraj , 'xLim', [0,1], 'xlab', '\theta',...
    'plotTitle', {['Prpsl.SD = ', num2str(proposalSD, '%g')], ...
    ['Eff.Sz. = ', num2str(neff,'%g')]});

%  Trajectory, a.k.a. trace plot, end of chain:
subplot(3,1,2);
idxToPlot = (trajLength-100):trajLength;
plot( trajectory(idxToPlot) , idxToPlot, '-o', ...
    'MarkerFaceColor',[0.4 0.7 1], 'MarkerSize',2);
title('End of Chain');
xlabel('\theta'); xlim([0,1]); ylabel('Step in Chain');
hold on;
% Display proposal SD and acceptance ratio in the plot.
mystr = sprintf('$$\\frac{N_{acc}}{N_{pro}}=%.2g$$',...
    nAccepted/length(acceptedTraj));
t = text( 0.1 , trajLength-25 , mystr);
set(t,'Interpreter','latex');
hold off;

%  Trajectory, a.k.a. trace plot, beginning of chain:
subplot(3,1,3);
 idxToPlot = 1:100;
plot( trajectory(idxToPlot) , idxToPlot, '-o', ...
    'MarkerFaceColor',[0.4 0.7 1], 'MarkerSize',2);
title('Beginning of Chain');
xlabel('\theta'); xlim([0,1]); ylabel('Step in Chain');

img = getframe(gcf);
imwrite(img.cdata,fullfile('figures','BernMetrop0.png'));
