% This script creates the field level confusion matrix in html format ('data/fieldlevel.html') and in latex format ('data/Table2.tex')

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

fprintf(1,'Creating the field level confusion matrix (for Table 2)...\n\n'); 

htmlstr = 'data/Table2.html';
fid = fopen(htmlstr, 'w');

for testclassNum=4:4
    strsave = sprintf('data/yeastfieldclassify_c%d',testclassNum);
    load([strsave '.mat']);
    rownames = cell(1,testclassNum);

    for t=1:testclassNum
        rownames{t} = [classNames{testorder(t)}];
    end

    fprintf(fid, 'Perfrom %d fold cross-validation on %d largest classes<BR>\n', nfold, testclassNum);
    
    fprintf(fid, 'Confusion matrix:');
    fprintf(fid, 'Accuracy(weighted by class):%.1f%%<BR>\n', 100*sum(diag(conf))/sum(conf(:)));
    confnorm = 100*conf./(sum(conf,2)*ones(1,testclassNum));
    fprintf(fid, 'Unweighted accuracy:%.1f%%<BR>\n', 100*sum(diag(confnorm))/sum(confnorm(:)));
   
    tablestr = tz_datatable(confnorm,rownames,rownames,'html');
    fprintf(fid, '%s<BR>\n', tablestr);
end

fclose(fid);


fprintf(1,'Creation of %s is done!\n\n', htmlstr);
