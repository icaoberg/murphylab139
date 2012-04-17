function display(ker)

% DISPLAY
%
% Display a textual representation of an exponential radial basis kernel object.
%
%    display(ker);

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

%
% File        : @exprbf/display.m
%
% Author      : Kai Huang
%
% Description : Display method for a class implementing an exponential radial basis kernel.
%               Part of an object-oriented implementation of Vapnik's Support
%               Vector Machine, as described in [1].  
%
% References  : [1] V.N. Vapnik,
%                   "The Nature of Statistical Learning Theory",
%                   Springer-Verlag, New York, ISBN 0-387-94559-8,
%                   1995.
%

fprintf(1, '\nexponential radial basis kernel (gamma = %f)\n\n', ker.gamma);

% bye bye...

