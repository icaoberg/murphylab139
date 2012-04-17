function ker = exprbf(arg)

% Exponential Radial Basis Function
%
% Construct an exponential radial basis kernel object,
%
%    K(x1, x2) = exp(-gamma.*(sum(abs(x1 - x2))))
%
% Examples:
%
%    % default constructor (spherical RBF kernel, gamma = 1.0)
%
%    ker1 = exprbf;
%
%    % copy constructor
%
%    ker2 = exprbf(ker1);
%
%    % construct exponential radial basis kernel, gamma = 0.5
%
%    ker3 = exprbf(0.5);

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
% File        : @exprbf/exprbf.m
%
% Description : Constructor for a class implementing an exponential
%               radial basis kernel,
%               forming part of a Matlab toolbox implementing Vapnik's
%               Support Vector Machine, as described in [1].
%
% References  : [1] V.N. Vapnik,
%                   "The Nature of Statistical Learning Theory",
%                   Springer-Verlag, New York, ISBN 0-387-94559-8,
%                   1995.
%

if nargin == 0
   
   % this is the default constructor
   
   ker.gamma = 1;
   ker       = class(ker, 'exprbf');
   
elseif isa(arg, 'exprbf');
   
   % this is the copy constructor
   
   net = arg;
   
elseif nargin == 1
   
   if prod(size(arg)) ~= 1

      error('gamma must be a scalar');

   end

   ker.gamma = 1/(2*(arg^2));
   ker       = class(ker, 'exprbf');

else

   % there are no other constructors

   help exprbf

end

% bye bye...




