% This script summarizes the classification result

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

load('data/yeastfeat_d1.mat');
unique_num = sum(yeastdata.ucsfresult.classSize);
gene_num = length(yeastdata.geneNames);
unvisualized_num = sum(yeastdata.ucsf.classLabels == 1);

length(unique(yeastdata.geneNames)) - unvisualized_num;   % 4155 ORFs
length(yeastdata.geneNames);   % 4239 images
visualized_num = gene_num - unvisualized_num;   % 4239 - unvisualized_num

ambiguous = load('data/yeastfeat_d41.mat');
ambiguous_num = ambiguous.yeastdata.ucsfresult.classSize;   % 156

mixture = load('data/yeastfeat_d2.mat');
mixture_num = sum(mixture.yeastdata.ucsfresult.classSize);  % 1294

unprocessed_num = gene_num - unvisualized_num - ambiguous_num - mixture_num - unique_num ;
processed_num = ambiguous_num + mixture_num + unique_num;

fprintf(1,'Among %d ORFs downloaded from the website, %d are unvisualized.\n', gene_num, unvisualized_num);
fprintf(1,'Among %d visualized ORFs, %d are not successfully processed.\n', visualized_num, unprocessed_num);
fprintf(1,'Among %d visualized ORFs, %d are successfully processed.\n', visualized_num, processed_num);
fprintf(1,'Among %d processed ORFs, %d belong to ambiguous class.\n',processed_num,ambiguous_num);
fprintf(1,'Among %d processed ORFs, %d has unique labels (other than ambiguous).\n',processed_num, unique_num);
fprintf(1,'Among %d processed ORFs, %d has combinations of labels.\n',processed_num, mixture_num);



