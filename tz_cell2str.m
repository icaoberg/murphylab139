function s=tz_cell2str(x,brk)
%TZ_CELL2STR Convert a cell array to a string.
%   TZ_CELL2STR(X) returns a string, which is the comma separated 
%   combination of strings or numbers in the cell array X.
%
%   TZ_CELL2STR(X,D) lets a user specify a customized sparator by D.
%
%   Example:
%       tz_cell2str({'apple','orange','banana'}) returns
%       'apple,orange,banana'
%       tz_cell2str({1,2,3}) returns '1,2,3'
%   
%       tz_cell2str({'apple','orange','banana'},'-') returns
%       'apple-orange-banana'
%
%   See also TZ_NUM2STR

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

%   23-MAY-2003 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin<1
    error('1 or 2 arguments are required.');
end

if ~iscell(x)
    error('The first argument must be a cell array.');
end

if isempty(x)
    s=[];
    return
end

if ~exist('brk','var')
    brk=',';
end

n=length(x);
s=[];
for i=1:n-1
    s=[s,num2str(x{i}),brk];
end
s=[s,num2str(x{end})];