function showclassify(flag, seedidx)
% SHOWCLASSIFY loads results of cell level classification and generates supplementary material 
%   for a list proteins in the ambiguous category (data/supplement2.txt)  (When flag == 41)
%   or a list proteins in the punctate_composite category (data/supplement3.txt) (When flag == 42)
%
%   FLAG specifies either ambiguous class or punctate_composite class
%   SEEDIDX specifies the number of seeds for plurality voting

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

fprintf(1,'Retrieving results of cell level classification...\n\n');
resultroot = [pwd filesep 'celllevel/result'];


yeastfeatstr = sprintf(['data' filesep 'yeastfeat_d%d.mat'],3);
tmp = load(yeastfeatstr);
uniquelabels = tmp.uniquelabels;
allclassLabels = tmp.allclassLabels;
strClassNames = tmp.strClassNames;
[imagenums,testorder]=sort(tmp.yeastdata.ucsfresult.classSize, 'descend');
clear tmp;

classnum = 20;
if flag == 41
    checkgenename = 'YFL034W';  % for Figure 5
    tmp = load([resultroot '/SEED_F41C20N6A1M140R' num2str(seedidx) '_S1W1_C11G5N0libsvm.mat']);   
    titlestr = 'List of computer-assigned labels for 156 proteins within the ambiguous category.';
    savefilestr = ['data' filesep 'supplement2.txt'];    
end

if flag == 42
    checkgenename = 'YGR130C';  % for Figure 4
    tmp = load([resultroot '/SEED_F42C20N6A1M140R' num2str(seedidx) '_S1W1_C11G5N0libsvm.mat']);   
    titlestr = 'List of computer-assigned labels for 72 proteins within the punctate_composite category.';
    savefilestr = ['data' filesep 'supplement3.txt'];
end
checkidx = strmatch(checkgenename, tmp.yeastdata.geneNames);

yeastdata = tmp.yeastdata;
votenum = seedidx;

nclassifier = votenum * ones(length(yeastdata.ucsf.locNums),1) * 6;
mcorrect = zeros(length(yeastdata.ucsf.locNums),1);
assignment = zeros(length(yeastdata.ucsf.locNums),1);

count = 0;
% Goal: assign an location integer to each class.  (first stage)
for i=1:1
    for j=1:yeastdata.ucsfresult.classSize(i)
        freq = zeros(1,classnum);
        for k=1:length(yeastdata.ucsfresult.imganot{i}{j}.voteclass)       
            freq(yeastdata.ucsfresult.imganot{i}{j}.voteclass(k)) = ...
               freq(yeastdata.ucsfresult.imganot{i}{j}.voteclass(k)) + 1;
        end
        [m, midx] = max(freq);

        % pularity      
        assignment(yeastdata.ucsfresult.classidx{i}(j)) = uniquelabels(testorder(midx));     
        mcorrect(yeastdata.ucsfresult.classidx{i}(j)) = m;

        if yeastdata.ucsfresult.classidx{i}(j) == checkidx
            freqclass =find(freq > 0);
            fprintf(1,'%s is classified into the following class:\n', checkgenename);
            for f=1:length(freqclass)
                fprintf(1,'%s  %.2f%%\n', strClassNames{uniquelabels(testorder(freqclass(f)))}, 100*freq(freqclass(f))/sum(freq));
            end                 
            fprintf(1,'\n');
        end
    end
end
confidence = mcorrect ./ nclassifier;
fid = fopen(savefilestr,'w');

fprintf(fid,'Supplementary material for Automated Image Analysis of Protein Localization in Budding Yeast\n');
fprintf(fid,'%s\n', titlestr);
fprintf(fid,'ORF\tAutomated localization\tConfidence\n'); 

processed_idx = [];
unprocessed_idx = [];


count = 0;
counthigh = 0;

[y, sortidx] = sort(confidence, 'descend');
for i=1:length(yeastdata.geneNames)
    if yeastdata.ucsf.classLabels(sortidx(i)) == 1
        % unvisulized one
        continue;
    end

    if allclassLabels(sortidx(i)) == assignment(sortidx(i))
         % The label agrees 
        continue;
    end

    if assignment(sortidx(i)) > 0
    % The label does not agree
          fprintf(fid,'%s\t', yeastdata.geneNames{sortidx(i)});
          fprintf(fid,'%s', strClassNames{assignment(sortidx(i))});
          processed_idx = [processed_idx; sortidx(i)];
          fprintf(fid,'\t%.2f%%', 100*confidence(sortidx(i)));
          fprintf(fid,'\n');
          count = count + 1;
          if confidence(sortidx(i)) == 1
              counthigh = counthigh + 1;
          end
     else
          unprocessed_idx = [unprocessed_idx; sortidx(i)];
     end
end
fclose(fid);

fprintf(1,'Creation of %s is done!\n', savefilestr);
fprintf(1,'%d of proteins (%d of them have 100%% confidence) are put in the list\n\n', count, counthigh);


