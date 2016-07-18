function[HDI] = HDIofICDF( ICDFname , credMass, pdf_params )
%Arguments:
%   ICDFname is R's name for the inverse cumulative density function
%     of the distribution.
%   credMass is the desired mass of the HDI region.
%   tol is passed to R's optimize function.
%   pdf_params is an array including parameters specific to the pdf
%   being used, such as mu, sigma, a, b, theta...
%Return value:
%   Highest density iterval (HDI) limits in a vector.
% Example of use: For determining HDI of a beta(30,12) distribution, type
%   HDIofICDF( qbeta , shape1 = 30 , shape2 = 12 )
%   Notice that the parameters of the ICDFname must be explicitly named;
%   e.g., HDIofICDF( qbeta , 30 , 12 ) does not work.
% Adapted and corrected from Greg Snow's TeachingDemos package.

if ~exist('credMass','var')
    credMass = 0.95;
end

incredMass =  1.0 - credMass;

if numel(pdf_params)==2
f = @(point) (icdf(ICDFname,point+credMass,pdf_params(1),pdf_params(2))-...
    icdf(ICDFname,point,pdf_params(1),pdf_params(2)));
elseif numel(pdf_params==1)
    f = @(point) (icdf(ICDFname,point+credMass,pdf_params)-...
    icdf(ICDFname,point,pdf_params));
else
    f = @(point) (icdf(ICDFname,point+credMass)-...
    icdf(ICDFname,point));
end

minpoint = fminbnd(f,0,incredMass);
HDI = [icdf(ICDFname,minpoint,pdf_params(1),pdf_params(2)) ...
    icdf(ICDFname,minpoint+credMass,pdf_params(1),pdf_params(2))];


end
% Kruschke, J. K. (2011). Doing Bayesian data analysis: A
% Tutorial with R and BUGS. Elsevier Science/Academic Press.