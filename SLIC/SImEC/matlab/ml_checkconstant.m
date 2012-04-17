function status=ml_checkconstant(s1,s2)
%ML_CHECKCONSTANT Check constant for vectors
%
%   STATUS=ML_CHECKCONSTANT(X) returns 0 if all numbers in X are the same. 
%   Otherwise it returns 1.
%
%   STATUS=ML_CHECKCONSTANT(X1,X2) returns the one of the following values:
%       0 - all the same
%       1 - both are contants, but x1 is greater than X2
%       -1 - both are constants, but X1 is less than X2
%       2 - X1 is constant, X2 is not
%       3 - X2 is constant, X1 is not
%       4 - neither one is constant

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

status=4;

mins1=min(s1);
mins2=min(s2);
maxs1=max(s1);
maxs2=max(s2);
if mins1==maxs1 & mins2==maxs2
    if mins1==mins2
        status=0;
    else
        status=sign(mins1-mins2);
    end
else
    if mins1==maxs1
        status=2;
    end
    if mins2==maxs2
        status=3;
    end
end
