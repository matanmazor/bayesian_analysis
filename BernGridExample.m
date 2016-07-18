Theta = 0:0.1:1;                    % Sparse teeth for Theta.
pTheta = min( Theta , 1-Theta );    % Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta);        % Make pTheta sum to 1.0
Data = [repelem(0,0),repelem(1,3)]; % Single flip with 1 head

figure('Units', 'Centimeters', 'Position', [1, 3, 10, 14]);
posterior = BernGrid( Theta, pTheta , Data);% , plotType='Bars' , 
                    %  showCentTend='one' , showHDI=0 , showpD=0 )
saveas(gcf, 'BernGridExample0.eps')