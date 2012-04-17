function tz_save(filepath,vars,varnames,script,comments,options)
%TZ_SAVE Save variables with related information
%   TZ_SAVE(FILEPATH,VARS,VARNAMES,SCRIPT,COMMENTS) save the variable in 
%   the cell  array VARS in the file FILEPATH. The names of the variables 
%   are in the [string array] VARNAMES. SCRIPT is the script which created
%   the variables. COMMENTS is a string of comments for the data.
%   
%   TZ_SAVE(FILEPATH,VARS,VARNAMES,SCRIPT,COMMENTS,OPTIONS) also specifies 
%   options allowed in the the function SAVE. OPTIONS is a [string array].
%   
%   See also

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

%   24-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 5
    error('5 or 6 arguments are required')
end

for i=1:length(varnames)
    eval([varnames{i} '= vars{i};']);
end
machine = tz_machineinfo('name','computer','domain');
savetime = clock;
matlab_version = version;

if ~exist('options','var')
    options = {};
end

save(filepath,varnames{:}, ...
    'script','comments','machine','savetime','matlab_version',...
    options{:});
