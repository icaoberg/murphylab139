function sc_classyeast_field(testclassNum, nfold)
%SC_CLASSYEAST_FIELD runs field level classification on unique location classes using nfold cross-validation scheme
%   TESTCLASSNUM specifies number of classes to be tested
%   NFOLD specifies number of folds for cross-validation
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

%   11-Apr-2007 Initial write S.C. Chen
%   Copyright (c) Center for Bioimage Informatics, CMU

whichdatabase = 2;
strsave = sprintf('data/yeastfieldclassify_c%d.mat',testclassNum);
load('data/yeastfieldfeat.mat');

% classNum
% classSize
% classNames
% featbyClass

[imgperclass,order]=sort(classSize);
order = fliplr(order);
imgperclass = fliplr(imgperclass);

trials = sum(imgperclass>=5);

conf = zeros(testclassNum,testclassNum);

testimgperclass = imgperclass(1:testclassNum);
testorder = order(1:testclassNum);

allfeatcell = cell(1,testclassNum);    
for t=1:testclassNum
    fprintf(1,'Loading %d class for SDA...\n',t);        
    for i=1:testimgperclass(t)
        %[t, i]
        allfeatcell{t} = [allfeatcell{t}; featbyclass{testorder(t)}{i}];
    end        
end
[idx_sda, ignoreidx] = ml_stepdisc(allfeatcell,['data' filesep 'field' num2str(whichdatabase) num2str(testclassNum) '.sdalog']);

conf = svm_cv(allfeatcell, idx_sda, nfold, '-s 0 -c 8192 -g 0.001953125');

confnorm = conf./(sum(conf,2)*ones(1,testclassNum));

fprintf(1,'accuracy (over the diagonal of the normalized confusion marix) = %.1f%%\n',100*sum(diag(confnorm))/sum(confnorm(:)));

confnorm
save(strsave);
