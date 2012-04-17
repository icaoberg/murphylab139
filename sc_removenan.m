function [newfeatmat, validrow, nanrow] = sc_removenan(featmat)
%SC_REMOVENAN remove NaN in a feature matrix
%   [NEWFEATMAT, VALIDROW, NANROW] = SC_REMOVENAN(FEATMAT) returns a feature matrix without NaN features.
%   VALIDROW is the index for valid rows in the original feature matrix.
%   NANROW is the index for rows with NaN features in the original feature matrix.
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

%   23-Jan-2006 Initial write A. Shanghag
%   11-Apr-2007 Modified S.C. Chen
%   Copyright (c) Center for Bioimage Informatics, CMU

[x y] = find(isnan(featmat) ~=0);
nanrow = unique(x);
validrow = setdiff(1:size(featmat,1) , nanrow);
newfeatmat = featmat(validrow', :);
