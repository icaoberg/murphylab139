function sv = getsv(net)

% GETSV
%
% Accessor method returning the support vectors of a support vector 
% classifier network.
%
%    sv = getsv(net);

%
% File        : @svc/getsv.m
%
% Date        : September 12th September 2000
%
% Author      : Dr Gavin C. Cawley
%
% Description : Part of an object-oriented implementation of Vapnik's Support
%               Vector Machine, as described in [1].
%
% References  : [1] V.N. Vapnik,
%                   "The Nature of Statistical Learning Theory",
%                   Springer-Verlag, New York, ISBN 0-387-94559-8,
%                   1995.
%
% History     : 07/07/2000 - v1.00
%               12/09/2000 - v1.01 minor improvements to comments and help
%                                  message
%
% Copyright   : (c) Dr Gavin C. Cawley, September 2000
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
%

sv = net.sv;

% bye bye...

