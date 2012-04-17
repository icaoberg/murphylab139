function setcvspath(cvsroot,cvsver)
% SETCVSPATH(cvsroot, cvsver)
% This function add all the paths of the SLIC using addpath command
%
% Input:
%   cvsroot: the root directory of the SLIC
%   cvsver:  the version number of SLIC; default value is 2
%
% Output:
%   None
%
%   Usage example:	setcvspath('SLIC');

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


if cvsroot(1) ~= '/'
    error('Please input full path.');
end

if ~exist('cvsver','var')
    cvsver=2;
end

switch cvsver
case 1
    addpath([cvsroot '/3D'])
    addpath([cvsroot '/3D/bin'])
    addpath([cvsroot '/3D/demo'])
    
    addpath([cvsroot '/3D/matlab'])
    
    addpath([cvsroot '/input/matlab'])
    addpath([cvsroot '/input'])
    addpath([cvsroot '/preprocess/matlab'])
    
    addpath([cvsroot '/segment/matlab'])
    addpath([cvsroot '/SImEC/matlab'])
    addpath([cvsroot '/SImEC/demo/matlab'])
    addpath([cvsroot '/classify/matlab'])
    addpath([cvsroot '/classify/svm_v0.54'])
    addpath([cvsroot '/classify/netlab/matlab'])
    addpath([cvsroot '/data'])
    addpath([cvsroot '/featcalc/matlab'])
    addpath([cvsroot '/featcalc/features'])
    addpath([cvsroot '/featcalc/demo'])
    addpath([cvsroot '/misc'])
    addpath([cvsroot '/SDA/matlab'])
    addpath([cvsroot '/shared/matlab'])
    
    addpath([cvsroot '/shared/bin'])
    addpath([cvsroot '/TypIC/matlab'])
    
    % addpath([cvsroot '/objclassif/matlab'])
    addpath([cvsroot '/TypIC/matlab/mex'])
    addpath([cvsroot '/shared/matlab/mex'])
    addpath([cvsroot '/segment/matlab/mex'])
    addpath([cvsroot '/3D/matlab/mex'])
    addpath([cvsroot '/input/matlab/mex'])
    
case 2
    addpath([cvsroot '/3D/matlab']);
    addpath([cvsroot '/3D/matlab/mex']);
    addpath([cvsroot '/classify/matlab']);
    addpath([cvsroot '/classify/svm_v0.54']);
    addpath([cvsroot '/classify/netlab/matlab']);
    if exist([cvsroot '/classify/libsvm'],'dir')
        addpath([cvsroot '/classify/libsvm']);
    end
    addpath([cvsroot '/clustering/matlab']);
    addpath([cvsroot '/featcalc/matlab']);
    addpath([cvsroot '/featcalc/matlab/mex']);
    addpath([cvsroot '/input/matlab']);
    addpath([cvsroot '/input/matlab/mex']);
    addpath([cvsroot '/misc']);
    addpath([cvsroot '/preprocess/matlab']);
    addpath([cvsroot '/SDA/matlab']);
    addpath([cvsroot '/segment/matlab']);
    addpath([cvsroot '/SImEC/matlab']);
    addpath([cvsroot '/TypIC/matlab']);
    addpath([cvsroot '/TypIC/matlab/mex']);
    addpath([cvsroot '/genmodel/matlab']);
    addpath([cvsroot '/genmodel/mixture']);
end
