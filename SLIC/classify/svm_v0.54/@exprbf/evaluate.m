function K = evaluate(ker, x1, x2)

% EVALUATE
%
% Evaluate an exponential radial basis kernel, for example
%
%    K = evaluate(kernel, x1, x2);
%
% where x1 and x2 are matrices containing input patterns, where each column
% represents a variable and each row represents an observation.

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
% File        : @exprbf/evaluate.m
%
% Author      : Kai Huang
%
% Description : Evaluate an exponential radial basis kernel.  Part of an object-oriented
%               implementation of Vapnik's Support Vector Machine, as
%               described in [1].  
%
% References  : [1] V.N. Vapnik,
%                   "The Nature of Statistical Learning Theory",
%                   Springer-Verlag, New York, ISBN 0-387-94559-8,
%                   1995.
%

n = size(x2,1);
m = size(x1,1);
K = zeros(m,n);

if (m <= n)
   for p = 1:m
      K(p,:) = sum(abs(ones(n,1)*x1(p,:) - x2),2).';
   end
else
   for p = 1:n
      K(:,p) = sum(abs(x1 - ones(m,1)*x2(p,:)),2);
   end
end

K = exp(-ker.gamma*K);

% bye bye...

