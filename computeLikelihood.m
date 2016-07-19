function pDataGivenTheta = computeLikelihood( theta , data)
%compute the Bernoulli likelihood for a vector of theta values 
pDataGivenTheta = (theta.^sum( data ) .* (1-theta).^sum( data==0));
pDataGivenTheta(theta > 1 | theta < 0 ) = 0;

end

