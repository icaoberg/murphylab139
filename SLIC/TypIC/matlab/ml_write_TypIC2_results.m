function ml_write_TypIC2_results (id, file, crop, featsetname, featselection, numpc, pcntexplained, robustcovar, indices, distances, imagenames)
%Taken from write_TypIX_results.m from Dr. Murphy but writes a web page for the results instead of a text file.
 
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

fid = fopen('results.html', 'w');
fprintf(fid, '<HTML><HEAD><TITLE>TypIC 2.0 Report</TITLE></HEAD>\n');
fprintf(fid, '<BODY aLink=#0077ff bgColor=#ffffff link=#118811 text=#000000 vLink=#ff0000>\n');  

fprintf(fid, '<CENTER><B><U>TypIC 2.0 report'); 
if (id ~= -1) 
     st_time = datestr(id/10000); 
     fprintf(fid, [' for TypIC 2.0 session started on&nbsp;' st_time ' GMT']);
end

fprintf(fid, '.</U></B></CENTER>\n<BR>\n\n'); 
fprintf(fid, '<P><CENTER>The <B>Typ</B>ical <B>I</B>mage <B>C</B>hooser has extracted features from the images and has ranked them in order of typicality.  For more information regarding TypIC 2.0 please click <A HREF="/projects/typical/">here</A>.</CENTER></P>\n');

fprintf(fid, '<DIV align=center><B><U>Parameter Information</U></B>\n');
fprintf(fid, '<TABLE align=center border=1>\n');
fprintf(fid, '<TR><TD><B>Protein Directory Name</B></TD>');
fprintf(fid, '%s%s%s\n', '<TD align = center>', file, '</TD></TR>');
if ((isempty(crop)) == 0)
     fprintf(fid, '<TR><TD><B>Crop Directory Name</B></TD>');
     fprintf(fid, '%s%s%s\n', '<TD align = center>', crop, '</TD></TR>');
end
fprintf(fid, '<TR><TD><B>Feature Set Used</B></TD>');
fprintf(fid, '%s%s%s\n', '<TD align = center>', featsetname, '</TD></TR>');
fprintf(fid, '<TR><TD><B>Feature Selection Method</B></TD>');
fprintf(fid, '%s%s%s\n', '<TD align = center>', featselection, '</TD></TR>');
if ((strcmp(featselection, 'princomp') ~= 0))
     fprintf(fid, '<TR><TD><B>Percent Variance Explained</B></TD>');
     fprintf(fid, '%s%5.2f%s\n', '<TD align = center>', pcntexplained, '</TD></TR>');
     fprintf(fid, '<TR><TD><B>Number of Principal Components Used</B></TD>');
     fprintf(fid, '%s%d%s\n', '<TD align = center>', numpc, '</TD></TR>');
end
fprintf(fid, '<TR><TD><B>Robust Covariance Used</B></TD>');
fprintf(fid, '%s%s%s\n', '<TD align = center>', robustcovar, '</TD></TR>');
fprintf(fid, '</TABLE></DIV>\n');

fprintf(fid, '<P><CENTER>Images are ranked on scale of typicality from 1 to 0 based on Mahalanobis distance results. Images are displayed in a descending order of typicality.</CENTER></P>\n');

num = size(indices);
maxfilelength = 0;
for i=1:num
  maxfilelength = max(maxfilelength,length(imagenames{i}));
end
fprintf(fid, '<HR>');
fprintf(fid, '<TABLE align=center border=1 cellpadding=2><TBODY>\n');
fprintf(fid, '<TR>\n');

fprintf(fid, '%s%s%s%s','<TD><B>RANK</B></TD>','<TD align=center><B>FILENAME</B></TD>','<TD><B>TYPICALITY</B></TD>','<TD><B>MAHALANOBIS DIST</B></TD>');
fprintf(fid, '</TR>\n');

mkdir('thumbs');
for i=1:num
  fprintf(fid, '<TR>\n');
  im = ml_readimage([file '/' imagenames{indices(i)}]);
  im = ml_preprocess(double(im), []);
  smallim = imresize(im, 0.5);
  
  imwrite(smallim, ['thumbs/' num2str(i) '.jpg'], 'jpg');
  fprintf(fid, '%s%d%s%s%s%d%s%s%s%s%10.7f%s%s%16e%s', '<TD align=right>', i, '</TD>' ,'<TD align=left>', ...
	       '<IMG SRC="./thumbs/', i, '.jpg"></IMG>', ...
	       imagenames{indices(i)}, '</TD>', '<TD align=right>', (num-i)/(num-1), '</TD>', ...
	       '<TD align=right>', distances(indices(i)), '</TD>');
  fprintf(fid, '</TR>\n');
end

fprintf(fid, '</TBODY></TABLE>');
fprintf(fid, '<BR><BR>');

fprintf(fid, '<P align=center>Thank you for using <b>TypIC 2.0</b>. </p><p align=center><b>TypIC 2.0</b> is a service of the Murphy Lab, Carnegie Mellon University, and is available at&nbsp <a href="http://murphylab.web.cmu.edu/services/TypIC2/">http://murphylab.web.cmu.edu/services/TypIC2</a>.</p>');
fprintf(fid, '</BODY></HTML>\n');
status = fclose(fid);

