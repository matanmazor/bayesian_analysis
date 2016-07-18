function[posterior] = BernGrid(Theta , pTheta , Data)

% Theta is vector of values between 0 and 1.
% pTheta is prior probability mass at each value of Theta
% Data is vector of 0's and 1's.

% Check for input errors:
if ( any( Theta > 1 | Theta < 0 ) )
    error('Theta values must be between 0 and 1')
end
if ( any( pTheta < 0 ) )
    error('pTheta values must be non-negative')
end
if abs(sum(pTheta)-1)>10^-9
    error('pTheta values must sum to 1.0')
end
if ~( all( Data == 1 | Data == 0 ) )
    error('Data values must be 0 or 1')
end

% Create summary values of Data
z = sum( Data ) % number of 1's in Data
N = length( Data )

% Compute the Bernoulli likelihood at each value of Theta:
pDataGivenTheta = Theta.^z .* (1-Theta).^(N-z)
%Compute the evidence and the posterior via Bayes' rule:
pData = sum( pDataGivenTheta .* pTheta )
pThetaGivenData = pDataGivenTheta.* pTheta / pData

%% Plot the results.

%% Plot the prior.
subplot(3,1,1)
yLim = [0,1.1*max(pTheta)];
bar( Theta, pTheta, 'EdgeColor','None','FaceColor',[0.4 0.7 1])
xlim([0,1]); ylim(yLim);
xlabel('\theta'); ylabel('p(\theta)'); title('prior');
meanTheta = sum( Theta .* pTheta ); hold on;
text(0.1, yLim(2)*0.8,['mean=',num2str(meanTheta)]);
modeTheta = Theta(find(pTheta== max( pTheta )));
text(0.1, yLim(2)*0.5,['mode=',num2str(modeTheta)]);
box off;
hold on;
credMass = 0.95;
HDIinfo = HDIofGrid( pTheta, credMass)
%plot HDI
HDI = [Theta(HDIinfo.indices(1)) Theta(HDIinfo.indices(end))];
line(HDI,[0,0],'Color','k','LineWidth',5); hold on;
text(mean(HDI),HDIinfo.height, sprintf('%d%% HDI',credMass*100),'HorizontalAlignment','center');
text(HDI(1),pTheta(HDIinfo.indices(1)), num2str(HDI(1)),'HorizontalAlignment','center');
text(HDI(2),pTheta(HDIinfo.indices(end)), num2str(HDI(2)),'HorizontalAlignment','center');

hold off;

%% Plot the likelihood: p(Data|Theta)

subplot(3,1,2)
yLim = [0,1.1*max(pDataGivenTheta)];
bar( Theta, pDataGivenTheta, 'EdgeColor','None','FaceColor',[0.4 0.7 1])
xlim([0,1]); ylim(yLim); 
xlabel('\theta'); ylabel('p(D|\theta)'); title('likelihood');
meanTheta = sum( Theta .* pDataGivenTheta/sum(pDataGivenTheta) ); hold on;
text(0.1, yLim(2)*0.7,['mean=',num2str(meanTheta)]);
modeTheta = Theta(find(pDataGivenTheta== max( pDataGivenTheta )));
text(0.1, yLim(2)*0.5,['mode=',num2str(modeTheta)]);
box off;
hold off;
                               
%% Plot the posterior.
subplot(3,1,3)
yLim = [0,1.1*max(pThetaGivenData)];
bar( Theta, pThetaGivenData, 'EdgeColor','None','FaceColor',[0.4 0.7 1])
xlim([0,1]); ylim(yLim);
xlabel('\theta'); ylabel('p(\theta|Data)'); title('posterior');
meanTheta = sum( Theta .* pThetaGivenData ); hold on;
text(0.1, yLim(2)*0.8,['mean=',num2str(meanTheta)]);
modeTheta = Theta(find(pThetaGivenData== max( pThetaGivenData )));
text(0.1, yLim(2)*0.5,['mode=',num2str(modeTheta)]);
box off;
hold on;
credMass = 0.95;
HDIinfo = HDIofGrid( pThetaGivenData, credMass)
%plot HDI
HDI = [Theta(HDIinfo.indices(1)) Theta(HDIinfo.indices(end))];
line(HDI,[0,0],'Color','k','LineWidth',5); hold on;
text(mean(HDI),HDIinfo.height, sprintf('%d%% HDI',credMass*100),'HorizontalAlignment','center');
text(HDI(1),pThetaGivenData(HDIinfo.indices(1)), num2str(HDI(1)),'HorizontalAlignment','center');
text(HDI(2),pThetaGivenData(HDIinfo.indices(end)), num2str(HDI(2)),'HorizontalAlignment','center');
posterior = pThetaGivenData;
hold off;
