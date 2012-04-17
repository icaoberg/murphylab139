function [train_norm, test_norm] = sc_featurenorm(train, test)
% SC_FEATURENORM - Normalize training and test data for LIBSVM traing/testing
%
% [TRAIN_NORM, TEST_NORM] = SC_FEATURENORM(TRAIN, TEST)
%    Normalizes the data in train to have range -1 to 1.
%    The values used to normalize train are also used to normalize test.
%
%   18-Jan-1999 Initial write M. Boland
%   11-Apr-2007 Modified S.C. Chen
%

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

% R - number of training samples
R = size(train,1) ;
% T - number of test samples
T = size(test,1) ;

train_norm = 2 * (train-(ones(R,1)*min(train))) ./ (ones(R,1)*(max(train)-min(train))) - 1;
if (T > 0)
  test_norm = 2 * (test-(ones(T,1)*min(train))) ./ (ones(T,1)*(max(train)-min(train))) - 1;
else
  test_norm = [] ;
end

