function ml_write_TypIX_results (indices,distances,imagenames)

% function ml_write_TypIX_results (indices,distances,imagenames)
% This function write the result of TypIC in a file named "result.txt"
% indices contains the rank of the images by typicality
% distances is the distances of images from the center in feature space
% imagenames are the names of the images

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

outfile=fopen('results.txt','w');
nfiles=length(indices);
maxfilelength=0;
for i=1:nfiles
 maxfilelength = max(maxfilelength,length(imagenames{i}));
end
oformat=strcat('%4s\t%',int2str(maxfilelength),'s\t%10s\t%16s\n');
fprintf(outfile,oformat,'RANK','FILENAME','TYPICALITY','MAHALANOBIS DIST');
fprintf(outfile,oformat,'====','========','==========','================');
oformat=strcat('%4i\t%',int2str(maxfilelength),'s\t%10.7f\t%16e\n');
for i=1:nfiles
 fprintf(outfile,oformat,i,imagenames{indices(i)},(nfiles-i)/(nfiles-1),distances(indices(i)));
end
fclose(outfile);
