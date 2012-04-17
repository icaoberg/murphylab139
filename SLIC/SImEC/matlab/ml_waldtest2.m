function [pvalue,ts] = ml_waldtest2(s1,s2)

%ML_WALDTEST2: the two-sample Wald test
%   PVALUE=ML_WALDTEST2 performs the Wald test on 
%   S1 and S2 to test if they have the same mean.
%   PVALUE is the p-value.
%   wald test: (X1bar-X2bar)/var ~ N(0,1)
%
%   Written by T. Zhao

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


status=ml_checkconstant(s1,s2);
switch status
case 0
    pvalue=1;
    ts=0;
case 1
    pvalue=0;
    ts=Inf;
case -1
    pvalue=0;
    ts=-Inf;
otherwise
    
    mean1=mean(s1);
    mean2=mean(s2);
    var1=var(s1,1);
    var2=var(s2,1);
    
    mean12=mean1-mean2;
    var12=var1/size(s1,2)+var2/size(s2,2);
    
    w=mean12/sqrt(var12);
    
    pvalue=2*normcdf(-abs(w),0,1);
    ts=w;
end