function [ok, measure]= sc_testmask(mask, dna)
%SC_TESTMASK tests if the mask is desired
%   [OK, MEASURE] = SC_TESTMASK(MASK, DNA)
%   
% Input:
%   MASK is the mask to be tested
%   DNA is the dna image
%
% Output:
%   OK is an flag which indicates if the mask pass the test
%   OK = 1  means it is a good mask
%   OK = -1 means the mask has bad eccentricity
%   OK = -2 means the mask is too small
%   OK = -3 means the mask is on the boundary of the image
%   MEASURE is a struct which stores the eccentricity, area and total dna
%
% usage:
%   [classok, measure] = sc_testmask(class==i, dnaimg);
%
%   See also: maintestmask
%
%   Nov-2005 initially written by G. Gordon
%   09-Aug-2006 Modified by S. C. Chen

mask = logical(round(mask+0));

% sizes
[n, m] = size(mask);
[xs, ys] = meshgrid(1:m, 1:n);

% get pixel coords and number of pixels
xm = xs(mask);
ym = ys(mask);
count = sum(mask(:));

mask = double(mask);

% mask with edges smoothly downweighted
bl = 11;
gaussfilt = -3:6/(bl-1):3;
gaussfilt = exp(-gaussfilt.^2/2);
gaussfilt = gaussfilt / sum(gaussfilt);
softmask = conv2(gaussfilt, gaussfilt, mask, 'same');
softmask = softmask .* mask;

% amount of DNA inside mask (counting edges less)
totdna = sum(sum(dna.*softmask));

% center of mass
mux = mean(xm);
muy = mean(ym);

% length of major/minor axes
sigma = svd([xm-mux, ym-muy], 0) ./ sqrt(count);

measure.sigma = min(sigma) / max(sigma);
measure.count = count;
measure.totdna = totdna;

ok = 1;

if (min(sigma) / max(sigma) < .5)
    fprintf('bad eccentricity %f\n', min(sigma) / max(sigma));
    ok = -1;
end

if (count < 249.8617 ) % mean - 2std
    fprintf('mask is too small %f\n', count);
    ok = -2;
end

if ((min(xm(:)) == 1) | (max(xm(:)) == m) | (min(ym(:)) == 1) | (max(ym(:)) == n))
    fprintf('the mask touchs the boundary of the image\n');    
    ok = -3;
end
