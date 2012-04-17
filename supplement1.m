% This script generates the confusion matrix for cell level classification (data/Table3.html)
% Draw precision-recall plot (data/Fig2.tif)
% ,and create Supplementary material for list of 501 proteins whose label from visual assignment differs from that by automated classification

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

  clear;


flag = 3;  seedidx = 25;
fprintf(1,'Retrieving results of cell level classification...\n\n');
resultroot = [pwd filesep 'celllevel/result'];

yeastfeatstr = sprintf(['data' filesep 'yeastfeat_d%d.mat'],flag);    
tmp = load(yeastfeatstr);
uniquelabels = tmp.uniquelabels;
allclassLabels = tmp.allclassLabels;
strClassNames = tmp.strClassNames;

clear tmp;
% flag 1 represents 21 class 
if flag == 1
    tmp = load([resultroot '/SEED_F1C21N6A1M140R' num2str(seedidx) '_S1W1_C11G5N0libsvm.mat']);
    classnum = 21;
end

% flag 3 represents 20 class (21 without punctate_composite class)
if flag == 3
    tmp = load([resultroot '/SEED_F3C20N6A1M140R' num2str(seedidx) '_S1W1_C11G5N0libsvm.mat']);
    classnum = 20; 
end

yeastdata = tmp.yeastdata;
testorder = tmp.testorder;

allidx = [];
for i=1:classnum
    allidx = [allidx; yeastdata.ucsfresult.classidx{i}];
end
allidx = sort(allidx, 'descend');

votenum = seedidx;
nclassifier = votenum * ones(length(yeastdata.ucsf.locNums),1);
mcorrect = zeros(length(yeastdata.ucsf.locNums),1);
assignment = zeros(length(yeastdata.ucsf.locNums),1);

YAL009Wclassidx = strmatch('YAL009W', yeastdata.geneNames);
% Goal: assign an location integer to each class.  (first stage)
for i=1:length(testorder)
    for j=1:yeastdata.ucsfresult.classSize(testorder(i))
        freq = zeros(1,classnum);
        for k=1:length(yeastdata.ucsfresult.imganot{testorder(i)}{j}.voteclass)       
            freq(yeastdata.ucsfresult.imganot{testorder(i)}{j}.voteclass(k)) = ...
               freq(yeastdata.ucsfresult.imganot{testorder(i)}{j}.voteclass(k)) + 1;
        end
        [m, midx] = max(freq);

        % pularity
        assignment(yeastdata.ucsfresult.classidx{testorder(i)}(j)) = uniquelabels(testorder(midx));
        mcorrect(yeastdata.ucsfresult.classidx{testorder(i)}(j)) = m;


        if yeastdata.ucsfresult.classidx{testorder(i)}(j) == YAL009Wclassidx 
    	    freqclass =find(freq > 0);
	    fprintf(1,'YAL009W is classified into the following class:\n');
            for f=1:length(freqclass)
                fprintf(1,'%s  %.2f%%\n', strClassNames{uniquelabels(testorder(freqclass(f)))}, 100*freq(freqclass(f))/sum(freq));            
            end
	    fprintf(1,'\n');
        end
    end
end
confidence = mcorrect ./ nclassifier;



savefilestr = ['data' filesep 'supplement1.txt'];    fid = fopen(savefilestr,'w');

fprintf(fid,'Supplementary material for Automated Image Analysis of Protein Localization in Budding Yeast\n');
fprintf(fid,'List of 501 proteins whose label from visual assignment differs from that by automated classification.\n');
fprintf(fid,'ORF\tUCSF localization\tAutomated localization\tConfidence\n'); 

processed_idx = [];
unprocessed_idx = [];

count = 0;
counthigh = 0;

correct = zeros(length(yeastdata.ucsf.locNums),1);

[y, sortidx] = sort(confidence, 'descend');
for i=1:length(sortidx)
    if yeastdata.ucsf.classLabels(sortidx(i)) == 1
        % unvisulized one
        continue;
    end

    if allclassLabels(sortidx(i)) == assignment(sortidx(i))
         % The label agrees 
	correct(sortidx(i)) = 1;
        continue;
    end

     if assignment(sortidx(i)) > 0
          % The label does not agree
          fprintf(fid,'%s\t%s', yeastdata.geneNames{sortidx(i)}, yeastdata.ucsf.strClassNames{yeastdata.ucsf.classLabels(sortidx(i))});
          fprintf(fid,'\t%s', strClassNames{assignment(sortidx(i))});
          processed_idx = [processed_idx; i];
          fprintf(fid,'\t%.2f%%', 100*confidence(sortidx(i)));
          fprintf(fid,'\n');
          count = count + 1;
	  if confidence(sortidx(i)) == 1          
              counthigh = counthigh + 1;
          end
     else
          unprocessed_idx = [unprocessed_idx; i];
     end
end
fclose(fid);

fprintf(1,'Creation of %s is done!\n', savefilestr);
fprintf(1,'%d of proteins (%d of them have 100%% confidence) are put in the list\n', count, counthigh);

conf = zeros(classnum, classnum);
for t=1:classnum
    for e=1:yeastdata.ucsfresult.classSize(testorder(t))
        c = find(uniquelabels(testorder) == assignment(yeastdata.ucsfresult.classidx{testorder(t)}(e)));
        conf(t,c) = conf(t,c) + 1;
    end
end

confnorm = conf./(sum(conf,2)*ones(1,classnum));
fprintf(1,'image accuracy = %.2f%% (average over the diagonal of confusion matrix: %.2f%%)\n',...
    100*sum(diag(conf))/sum(conf(:)), 100*sum(diag(confnorm))/sum(confnorm(:))  );

fprintf(1,'%d (%.2f%%) images are classified correct\n\n', sum(correct), 100*sum(correct)/sum(tmp.yeastdata.ucsfresult.classSize));


label = correct(allidx);
cellconfidence = confidence(allidx);
[precision, recall] = precision_recall(label, cellconfidence);

fid = fopen(['data' filesep 'Table3.html'], 'w');
for t=1:classnum
    rownames{t} = [yeastdata.ucsfresult.classNames{testorder(t)}];
end
tablestr = tz_datatable(confnorm,rownames,rownames,'html');
fprintf(fid, '%s', tablestr);
fclose(fid);
fprintf(1,'Confusion matrix file %s is created\n\n',['data' filesep 'Table3.html']);
