function[posterior]= BernBeta( priorBetaAB , Data)

  % priorBetaAB is two-element vector of beta a,b shape parameters
  % Data is vector of 0's and 1's.
  % source("DBDA2E-utilities.R") # for HDIofICDF()
  % Check for input errors:
  if  any( priorBetaAB < 0 )  
    error('priorBetaAB values must be nonnegative')
  end
  if length( priorBetaAB ) ~= 2  
    error('priorBetaAB must be a vector of two values')
  end
  if ~all( Data == 1 | Data == 0 ) 
    error('Data values must be 0 or 1')
  end
  
  % For notational convenience, rename components of priorBetaAB:
  a = priorBetaAB(1);
  b = priorBetaAB(2);
  
  % Create summary values of Data:
  z = sum( Data ); % number of 1's in Data
  N = length( Data ); 
  
  Theta = 0.001:0.001:0.999;   % points for plotting
  pTheta = betapdf( Theta , a , b );% prior for plotting
  pThetaGivenData = betapdf( Theta , a+z , b+N-z ); %posterior for plotting
  pDataGivenTheta = Theta.^z .* (1-Theta).^(N-z); % likelihood for plotting

  %% Plot the results.

%% Plot the prior.
subplot(3,1,1)
yLim = [0,1.1*max(pTheta)];
bar( Theta, pTheta, 'EdgeColor','None','FaceColor',[0.4 0.7 1])
xlim([0,1]); ylim(yLim);
xlabel('\theta'); ylabel(['betapdf(\theta|',num2str(a),',',num2str(b),')']); 
title('Prior (beta)');
meanTheta = sum( Theta .* pTheta /sum(pTheta)); hold on;
text(0.1, yLim(2)*0.8,['mean=',num2str(meanTheta)]);
modeTheta = Theta(find(pTheta== max( pTheta )));
text(0.1, yLim(2)*0.5,['mode=',num2str(modeTheta)]);
box off;
hold on;
credMass = 0.95;
HDI = HDIofICDF( 'beta',0.95,priorBetaAB);
% %plot HDI
HDI_height = betapdf(HDI(1),a,b);
plot(HDI,[HDI_height,HDI_height],'Color','k','LineWidth',0.5); hold on;
plot([HDI(1),HDI(1)],[0,HDI_height],':','Color','k'); hold on;
plot([HDI(2),HDI(2)],[0,HDI_height],':','Color','k'); hold on;
text(mean(HDI),HDI_height*1.7, sprintf('%d%% HDI',credMass*100),'HorizontalAlignment','center');
text(HDI(1),HDI_height*1.5, num2str(HDI(1),2),'HorizontalAlignment','center');
text(HDI(2),HDI_height*1.5, num2str(HDI(2),2),'HorizontalAlignment','center');

hold off;

%% Plot the likelihood: p(Data|Theta)

subplot(3,1,2)
yLim = [0,1.1*max(pDataGivenTheta)];
bar( Theta, pDataGivenTheta, 'EdgeColor','None','FaceColor',[0.4 0.7 1])
xlim([0,1]); ylim(yLim); 
xlabel('\theta'); ylabel('p(D|\theta)'); title('likelihood (Bernoulli)');
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
xlabel('\theta'); ylabel(['betapdf(\theta|',num2str(a+z),',',num2str(b+N-z),')']);
title('Posterior (Beta)');
meanTheta = sum( Theta .* pThetaGivenData/sum(pThetaGivenData) ); hold on;
text(0.1, yLim(2)*0.8,['mean=',num2str(meanTheta)]);
modeTheta = Theta(find(pThetaGivenData== max( pThetaGivenData )));
text(0.1, yLim(2)*0.5,['mode=',num2str(modeTheta)]);
box off;
hold on;
credMass = 0.95;
HDI = HDIofICDF( 'beta',0.95,[a+z , b+N-z]);
% %plot HDI
HDI_height = betapdf(HDI(1),a+z , b+N-z);
plot(HDI,[HDI_height,HDI_height],'Color','k','LineWidth',0.5); hold on;
plot([HDI(1),HDI(1)],[0,HDI_height],':','Color','k'); hold on;
plot([HDI(2),HDI(2)],[0,HDI_height],':','Color','k'); hold on;
text(mean(HDI),HDI_height*1.7, sprintf('%d%% HDI',credMass*100),'HorizontalAlignment','center');
text(HDI(1),HDI_height*1.5, num2str(HDI(1),2),'HorizontalAlignment','center');
text(HDI(2),HDI_height*1.5, num2str(HDI(2),2),'HorizontalAlignment','center');
posterior = pThetaGivenData;
hold off;

end