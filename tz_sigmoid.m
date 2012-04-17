function y = tz_sigmoid(x,threshold,slope)
%TZ_SIGMOID Sigmoid function.
%   Y = TZ_SIGMOID(X) returns 1/(1+e^(-X)).
%   
%   Y = TZ_SIGMOID(X,THRESHOLD) returns 1/(1+e^(threshold-X)).

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

%   09-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1-3 arguments are required')
end

if nargin < 2
    threshold = 0;
end

if nargin < 3
    slope = 1;
end

y = 1./(1+exp(slope*(threshold-x)));