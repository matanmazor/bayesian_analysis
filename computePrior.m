function [ pTheta ] = computePrior( theta )
%Define the prior density function. 
  pTheta = betapdf( theta , 1 , 1 );
  %The theta values passed into this function are generated at random,
  %and therefore might be inadvertently greater than 1 or less than 0.
  %The prior for theta > 1 or for theta < 0 is zero:
  pTheta(theta > 1 | theta < 0) = 0;
end



