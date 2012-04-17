function confusion_matrix = svm_cv(features, idx_sda, nfold, options)
%SVM_CV calculates nfold cross-validation with selected SDA index using LIBSVM
%   CONFUSION_MATRIX = SVM_CV(FEATURES, IDX_SDA, NFOLD, OPTIONS)
%   FEATURES is a cell-array which stores the matrix
%   IDX_SDA is the pre-selected SDA index
%   NFOLD specifies the number of fold for cross-validation 
%   OPTIONS specifies the parameters for LIBSVM
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

%   14-Mar-2007 Written by S.-C. Chen
%   Copyright (c) Center for Bioimage Informatics, CMU

k = size(features,2);
imagenums = zeros(k,1);

for i=1:k
    features{i} = double(features{i}(:,idx_sda));
    features{i} = [features{i} i*ones(size(features{i},1),1)];
    imagenums(i) = size(features{i},1);
end

confusion_matrix = zeros(k);
correct = 0;
rand('state', 0);

for j=1:k
    xo{j}  = features{j};
    idx{j} = randperm(imagenums(j));
    binsize = floor(imagenums(j)/nfold);

    start_idx = 1;
    for w=1:(nfold-mod(imagenums(j), nfold))
        data{j}{w} = xo{j}(idx{j}(start_idx:(start_idx+binsize-1)),1:end);
        start_idx  = start_idx+binsize;
    end

    for w=(nfold+1-mod(imagenums(j), nfold)):nfold
        data{j}{w} = xo{j}(idx{j}(start_idx:(start_idx+binsize)),1:end);
        start_idx  = start_idx+binsize+1;
    end
end

for w=1:nfold % n-fold cross validations

    x  = [];
    xt = [];
    yo = [];
    yot = [];

    for j=1:k

        sym = ones(1,nfold);
        sym(w) = 0;
        sym_idx = find(sym);

        for q = 1:length(sym_idx)
            x   = [ x; data{j}{sym_idx(q)}(:, 1:(end-1))];
            yo  = [yo; data{j}{sym_idx(q)}(:, end)];
        end

        xt = [xt; data{j}{w}(:,1:(end-1))];
        yot = [yot; data{j}{w}(:,end)];

        % the data matrix contains n columns where n is the number of input variable
        % each row represents a pattern and each column a variable.
        % the target matrix contains k columns where k is the number of classes
        % each row represents a pattern and each column a class. -1 means false
        % +1 means true.

    end

    [x xt] = sc_featurenorm(x, xt);
    model = svmtrain(yo, x, options);  
    [O, accuracy, dec_values] = svmpredict(yot, xt, model);
    for i=1:k
        for j=1:k
            confusion_matrix(i,j) = confusion_matrix(i,j)+length(find(yot == i & O == j));
        end
    end

end

for i=1:k
    for j=1:k
        if (j == i)
            correct = correct + confusion_matrix(j, j);
        end
    end
end

confusion_matrix
fprintf(1, 'correct=%d, the accuracy is %4.1f%%\n', correct, ...
    100*correct/sum(confusion_matrix(:)));

fprintf(1,'bye bye...\n');
