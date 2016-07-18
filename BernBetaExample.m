% Specify the prior:
t = 0.65;             % Specify the prior MODE.
n = 25;               % Specify the effective prior sample size.
a = t*(n-2) + 1;      % Convert to beta shape parameter a.
b = (1-t)*(n-2) + 1;  % Convert to beta shape parameter b.

Prior = [a,b];       % Specify Prior as vector with the two shape parameters.

% Specify the data:
N = 330;                         % The total number of flips.
z = 130;                         % The number of heads.
% Convert N and z into vector of 0's and 1's.
Data = [repelem(0,N-z),repelem(1,z)];
figure('Units', 'Centimeters', 'Position', [1, 3, 10, 14]);
posterior = BernBeta( Prior, Data);
saveas(gcf, 'BernBetaExample.png');
