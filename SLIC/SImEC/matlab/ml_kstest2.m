function [pvalue,ts]=ml_kstest2(s1,s2)

%ML_KSTEST2: Two-sample Kolmogorov-Smirnov goodness-of-fit hypothesis test.
%   [PVALUE,TS]=ML_KSTEST2 performs a Kolmogorov-Smirnov (K-S) test 
%   to determine if independent random samples, X1 and X2, are drawn from 
%   the same underlying continuous population.
%   The funciton will return p-value PVALUE and test statistic TS. 
%   It is the same with matlab function kstest2 except the return values.
%   See also
%   kstest2

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

[h,pvalue,ts]=kstest2(s1,s2);