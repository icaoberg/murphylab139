function index=ml_rank(A)

%ML_RANK: rank each row
%   INDEX=ML_RANK(A) gives rank for every element
%   in matirx A along rows and returns a rank matrix
%   INDEX.
%   For example, if A is [4 1 2;3 4 2;4 7 1], then
%   INDEX is [3 1 2;2 3 1;2 3 1]
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

[m,n]=size(A);
[B,idx]=sort(A,2);
for i=1:m
    index(i,idx(i,:))=1:n;
end