function feat = sc_maskfeatcalc(maski, dnaim)
%SC_MASKFEATCALC calculates features of a mask for post processing
%   FEAT = SC_MASKFEATCALC(MASKI, DNAIM)
%   
% Input:
%   MASKI is one mask of a single cell region (logical array)
%   DNAIM is the DNA image (double)
%
% Output:
%   FEAT is a 1 X 10 double array which stores the 10 features related to the mask
%
% usage: 
%   feat = sc_maskfeatcalc(maski, dnaim);
%
%   See also: sc_cutcells

%   09-Aug-2006 Modified by S. C. Chen

feat = zeros(1,10);
[xc, yc] = find(maski==1);
x = [xc yc];	meanx = mean(x);
A = (x - repmat(meanx, [size(x,1) 1]))'*(x - repmat(meanx, [size(x,1) 1]));
dist = sqrt(diag((x - repmat(meanx, [size(x,1) 1]))*inv(A)*(x - repmat(meanx, [size(x,1) 1]))'));
[eigvec,eigval] = eig(A);
eigval = diag(eigval);
[eigval, idx] = sort(eigval);	%sort in descending order
a = sqrt(eigval(2));  b = sqrt(eigval(1));

feattotpixel = sum(maski(:));
featellipsearea = pi * a * b * max(dist)^2;
featareafrac = feattotpixel / featellipsearea;
feateccentricity = sqrt(1-(b/a)^2);

feat(1) = feateccentricity;         % eccentricity: The eccentricity of the mask
feat(2) = featareafrac;             % area fraction: The area ratio of the mask and fitted ellipse
feat(3) = feattotpixel;             % total pixel: Total number of pixels
feat(5) = sum(sum(dnaim.*maski));   % DNA per pixel: Average DNA intensity per pixel in the mask
feat(4) = feat(5)/feat(3);          % total DNA: Total number of DNA intensity
feat(6:10) = feat(1:5).^2;          % synthesized: Square of each feature above
