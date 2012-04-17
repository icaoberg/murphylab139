function [precision, recall] = precision_recall(label, confidence)
% PRECISION_RECALL draws the precision recall plot
%   [PRECISION, RECALL] = PRECISION_RECALL(LABEL, CONFIDENCE)  
%
% Input:
%   LABEL specify the positive example 1 or negative example 0
%   CONFIDENCE specify the percentage of classifiers agrees
%
% Output:
%   PRECISION is a vector of precision value at different thresholds
%   RECALL is a vector of recall value at different thresholds
%   data/Fig2.tif
%
%   See also supplement1
%
%   06-June-2007 Written by S. C. Chen


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


[sortconfidence,sortidx] = sort(confidence,'descend');
sortlabel = label(sortidx);

totalpositive = sum(label);

tp=0;   fp=0;
fmeasure = [];
precision=[];
recall=[];

thres = sort(unique(confidence),'descend');

for i=1:length(thres)
    num = sum(confidence >= thres(i));   
    tp = sum(sortlabel(1:num)==1);
    fp = sum(sortlabel(1:num)==0);

    pn = tp/(tp+fp);
    rn = tp/totalpositive;
    precision(i) = pn;
    recall(i) = rn;
    fmeasure(i) = 2*pn*rn/(pn+rn);
end
 

[maxf, fidx] = max(fmeasure);

fprintf(1,'\nCalculating precision and recall...\n');
[thres';precision;recall]

h=figure(1);
plot(recall, precision, 'ko-','MarkerSize',5,'MarkerFaceColor',[0 0 0], 'LineWidth',2.5);
xlabel('Recall','FontSize',15);
ylabel('Precision','FontSize',15);
axis([0.7 1 0.7 1]);
p=get(h);
set(p.Children,'FontSize',15);
set(p.Children,'LineWidth',2);
print('-dtiff','-r300',['data' filesep 'Fig2']);
