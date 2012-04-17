function [F,df1,df2]=ml_fvalue(X,Y)

%ML_FVALUE: F-value between two matrices
%   [F,DF1,DF2]=ML_FVALUE(X,Y) calculates F-value between X and Y.
%   F is the F-value and DF1,DF2 are the degrees of the first and 
%   second freedom repectively.
%
%   See also
%   tz_ht2test2
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

p=size(X,2);
m=size(X,1);
n=size(Y,1);

%Calculate the Mahalanobsis distance between X and Y
dist=ml_twomaha(X,Y,1,1);

F=dist*(m+n-p-1)/(m+n-2)/p;
df1=p;
df2=m+n-p-1;