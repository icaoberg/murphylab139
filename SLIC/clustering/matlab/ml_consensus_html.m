function g = ml_consensus_html(feat, name, urls, htmlfile, html_URL, left,...
			       width, distance, method, penc, niter)

% FUNCTION G = ML_CONSENSUS_HTML(FEAT, NAME, URLS, HTMLFILE, HTML_URL, LEFT,...
%				 WIDTH, DISTANCE, METHOD, PENC, NITER)
% ml_consensus_html will perform k-means algorithm (if not done) followed by 
% AIC or BIC selection to choose the optimal k.  Then a set of dendrograms 
% built on randomly divided half set were created and used as input for
% majority consensus analysis. Finally a figure showing the consensus tree will
% be displayed and a corresponding web interface is created (provided the
% folder where htmlfile is created is configured to run cgi script
% feat: the input feature values. It is a cell array where each element is
%       the feature matrix for one clone
% name: name of each clone, a cell array with strings
% urls: urls for each image in the dataset.  A cell array, each element
%       representing a clone, each lement itself is a cell array where elements
%       are urls for individual images.  Same order as feat
% htmlfile: filename to store the final html output.  Should be ended with 
%           '.html'
% html_URL: corresponding url for htmlfile.
% left, width: parameters to calculate the position of the box corresponding
%              each leaf.  The fomula are:
%              Center of the box: left + idx_clone * width / total_number_clone
%              Width of the box: width/total_number_clone
% distance: specify the distance function.  One of the following values: 
%           'euclidean' (by default), 'mahalanobis'.
% method: Either 'aic' or 'bic'.  Default 'aic'.
% penc: pencentage of the threshold.  If a cluster has the highest number of
%       observations from a clone and the pencentage is greater than the penc,
%       The clone will be considered to be in this cluster. 33 by default.
% niter: Number of iterations in creating half set dendrograms (100 by default,
%        which yields 200 dedrograms).
% g: the resturn membership array from k-means/AIC(BIC).  Each row presents a
%    clone and first k columns stands for the number of the cells in that
%    cluster.  The k+1 column is the highest percentage of cells for that clone
%    in any of the clusters.  The k+2 column is the ID of the cluster which has
%    that largest percentage, if it is not less than perc. 
% Xiang Chen, Jan 04, 2005

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

  if (~exist('feat', 'var'))
     error('ml_consensus_html must have the argument: feat.  Type "help ml_consensus_html" for more information');
  end
  if (~exist('name', 'var'))
     for m = 1 : length(feat);
         name{m} = num2str(m);
     end
  end
  if(~exist('htmlfile', 'var'))
     htmlfile = ['consensustmp' num2str(length(feat)) '.html'];
  end
  currentdir = pwd;
  if (htmlfile(1) ~= '/')
	htmlfile = [currentdir '/' groupfile];
  end
  dots = findstr(htmlfile, '.');
  if (isempty(dots) | ~(strcmp(htmlfile(dots(end):end), '.html') | ...
			strcmp(htmlfile(dots(end):end), '.htm')))
     prefix = htmlfile;
     htmlfile = [htmlfile '.html'];
  else
     prefix = htmlfile(1:(dots(end) - 1));
  end
  dots = findstr(html_URL, '.');
  if (isempty(dots) | ~(strcmp(html_URL(dots(end):end), '.html') | ...
			strcmp(html_URL(dots(end):end), '.htm')))
     htmlprefix = html_URL;
     html_URL = [html_URL '.html'];
  else
     htmlprefix = html_URL(1:(dots(end) - 1));
  end
  groupfile = [prefix '.group.mat'];
  if (~exist('niter', 'var')|~isnumeric(niter))
     niter = 100;
  end
  if (~exist('distance', 'var'))
     distance = 'eu';
  end

  if (exist('method', 'var') & strcmp(method, 'bic'))
     bic = 1;
  else
     bic = 0;
     method = 'aic';
  end

  switch distance(1:2)
   case 'eu'
    eu = 1;
    ma = 0;
   case 'ma'
    ma = 1;
    eu = 0;
   otherwise
    error ('Method must be one of the following ''euclidean'', ''mahlanobis'' or ''both''.  Type ''help ml_consensus'' for details.');
end

  if (~exist('penc') | ~ isnumeric(penc))
     penc = 33;
  end

  if (~exist(groupfile, 'file'))
      ml_kmeans_aicbic(feat, groupfile, distance);
  end
  [g, c] = ml_readaicbic(groupfile, distance, method, penc);
  features = [];
  for m = 1 : length(feat)
          if (c{m})
             features{length(features)+1} = double(feat{m}(c{m}, :));
	     gene{length(features)} = name{m};
	     url{length(features)} = urls{m}(c{m});
	     %c2{length(features)} = c{m};
	  end
  end

  %if (~exist([prefix '.den'], 'file'))
      ml_bootstrap_half([prefix '.den'], '', features, gene, niter, ma);
      %end
 % if (~exist([prefix '.den_MRC'], 'file'))
      javadir = which('ml_consensus');
      javadir = [javadir(1:end-21) 'bin'];
      cd(javadir);
      unix(['java MajConsensus ' prefix '.den']);
      cd(currentdir);
      %end
  [links, outarr] = ml_readconsensus([prefix '.den_MRC'], features, ma);
  eval('label = ml_drawconsensus(links, gene, 9, ma);', ...
       'label = ml_drawconsensus(links(1:end-1), gene, 9, ma);');
  saveas(gcf, [prefix '.jpg'], 'jpeg');

  thresh = 0; %max number of images displayed
  init = 5; %init number of images displayed

     %htmlfile = ['/home/xiangc/cgi-bin/' img(1:end -3) 'html'];
     fid = fopen(htmlfile, 'w');
     if (fid == -1)
         error(['Failure in creating ' htmlfile '!']);
     end

     text = '<HTML>\n    <HEAD>\n        <TITLE>CONSENSUS TREE</TITLE>\n    </HEAD>\n    <BODY>\n        <P>Please click on the leaves of the consensus tree to see the representative images of that clone.\n        <BR><BR>\n        <P><IMG src="';
     text = [text htmlprefix '.jpg" width=1024 height=570 usemap="#scroll">\n        <MAP name="scroll">\n'];
     for m = 1 : length(label)
	     p1 = left + round((m - 0.5) * width / length(label));
	     p2 = left + round((m + 0.5) * width / length(label));
	     text = [text '            <area href="' htmlprefix '_clone.cgi?id=' num2str(m - 1) '&quan=' num2str(init) '" target="_blank" shape="rect"\n                  coords="' num2str(p1) ', 0, ' num2str(p2), ', 569">\n'];
     end

     if(~thresh)
         text = [text '        </MAP>\n        <form action="' htmlprefix '_quan.cgi" method="GET"\n            <P>Number of images shown for each clone: <input type="text" name="quan" value="5" size=2><br>\n            <input type="submit" value="Submit"><P>\n        </form>\n    </BODY>\n</HTML>\n'];
     else
         text = [text '        </MAP>\n        <form action="' htmlprefix '_quan.cgi" method="GET"\n            <P>Number of images shown for each clone (Maximum ' num2str(thresh) '): <input type="text" name="quan" value="' num2str(init) '" size=2><br>\n            <input type="submit" value="Submit"><P>\n        </form>\n    </BODY>\n</HTML>\n'];
     end
     fprintf(fid, text);
     fclose(fid);

     cgifile = [prefix '_quan.cgi'];
     fid = fopen(cgifile, 'w');
     if (fid == -1)
         error(['Failure in creating ' cgifile '!']);
     end
     
     if (~thresh)
         text = '#!/usr/bin/perl -wT\nuse CGI qw(:standard);\nuse strict;\n\nprint header;\nprint start_html("Consensus Trees");\n\n my $quan = param(''quan''); \n';
     else
         text = ['#!/usr/bin/perl -wT\nuse CGI qw(:standard);\nuse strict;\n\nprint header;\nprint start_html("Consensus Trees");\n\n my $quan = param(''quan''); \n if ($quan > ' num2str(thresh) ') {\n     $quan = ' num2str(thresh) ';}\n'];
     end
     text = [text 'print "        <P>Please click on the leaves of the consensus tree to see the representative images of that clone.\n        <BR><BR>\n        <P><IMG src=\\"' htmlprefix  '.jpg\\" width=1024 height=570 usemap=\\"#scroll\\">\n        <MAP name=\\"scroll\\">\n'];
     for m = 1 : length(label)
	     p1 = left + round((m - 0.5) * width / length(label));
	     p2 = left + round((m + 0.5) * width / length(label));
	     text = [text '            <area href=\\"' htmlprefix '_clone.cgi?id=' num2str(m - 1) '&quan=$quan\\" target=\\"_blank\\" shape=\\"rect\\"\n                  coords=\\"' num2str(p1) ', 0, ' num2str(p2), ', 569\\">\n'];
     end

     if(~thresh)
         text = [text '        </MAP>\n        <form action=\\"' htmlprefix '_quan.cgi\\" method=\\"GET\\"\n            <P>Number of images shown for each clone: <input type=\\"text\\" name=\\"quan\\" value=\\"$quan\\" size=2><br>\n            <input type=\\"submit\\" value=\\"Submit\\"><P>\n        </form>\n    </BODY>\n</HTML>\n"'];
     else
         text = [text '        </MAP>\n        <form action=\\"' htmlprefix '_quan.cgi\\" method=\\"GET\\"\n            <P>Number of images shown for each clone (Maximum ' num2str(thresh) '): <input type=\\"text\\" name=\\"quan\\" value=\\"$quan\\" size=2><br>\n            <input type=\\"submit\\" value=\\"Submit\\"><P>\n        </form>\n    </BODY>\n</HTML>\n"'];
     end
     fprintf(fid, text);
     fclose(fid);
     unix(['chmod 755 ' cgifile]);

     cgifile = [prefix '_clone.cgi'];
     fid = fopen(cgifile, 'w');
     if (fid == -1)
         error(['Failure in creating ' cgifile '!']);
     end
     
     text = '#!/usr/bin/perl -wT\nuse CGI qw(:standard);\nuse strict;\n\nprint header;\nmy %%form;\nmy @genes = (';
     text = [text '"' label{1} '"'];
     for m = 2 : length(label) 
	 text = [text ', "' label{m} '"'];
     end
     text = [text ');\n\n'];
     text = [text 'my $id = param("id");\n'];
     text = [text 'my $gene = $genes[$id];\n'];
     text = [text 'print start_html("Representative images for $gene");\n\n'];
     fprintf(fid, text);
    
     feat_tmp = [];
     for m = 1 : length(features)
         feat_tmp = [feat_tmp; features{m}];
     end

     stdf = std(feat_tmp);
     stdf(find(stdf==0)) = 1;
     if (ma)
         feat_tmp = feat_tmp ./ repmat(stdf, [size(stdf, 1) 1]);
         inv_cov = inv(cov(feat_tmp));
     end


     text = 'my %%nimg = (#';
     for m = 1 : length(features)
         text = [text ',\n\t\t"' gene{m} '", "' num2str(length(c{m})) '"'];
     end
     text = [text ' );\n'];
     fprintf(fid, text);

     text = 'my %%img = (#';
     for m = 1 : length(features) 
	 feat_tmp = features{m} ./ repmat(stdf, [size(features{m}, 1) 1]);
	 mf = mean(feat_tmp);
	 dist = [];
	 for n = 1 : size(feat_tmp, 1)
             if (ma)
                 dist(n) = (feat_tmp(n, :) - mf) * inv_cov * (feat_tmp(n, :) - mf)';
	     else
		 dist(n) = sum(((feat_tmp(n, :) - mf)).^2);
	     end
	 end
	 [V, I] = sort(dist);
	 for n = 1 : length(I)
	     text = [text ',\n\t\t"' gene{m} num2str(n) '", "' url{m}{I(n)} '"'];
	 end
     end
     text = [text ' );\n'];
     fprintf(fid, text);

     text = 'my $quan = param("quan");\n';
     fprintf(fid, text);

     text = 'if ($quan > $nimg{$gene}) {\n';
     text = [text '\t$quan = $nimg{$gene};\n}\n'];
     text = [text 'my @imgs = ();\n'];
     text = [text 'for (my $i = 1; $i <= $quan; $i++) {\n'];
     text = [text '\tmy $imgno = int(($i - 1) * $nimg{$gene} / $quan + 1);\n'];
     text = [text '\tmy $imgname = "$gene$imgno";\n'];
     text = [text '\tpush(@imgs, $img{$imgname});\n}\n'];
     fprintf(fid, text);

     text = 'print ''<CENTER>'';\n';
     text = [text 'print ''<TABLE border="0" width="900" cellspacing="0" cellpadding="0">'';\n'];
     text = [text 'print ''  <TBODY>'';\n'];
     text = [text 'print ''    <TR>'';\n'];
     text = [text 'print ''      <TD bgcolor="#999999" align="center" height="38"><p> <h2>Images</h2></TD>'';\n'];
     text = [text 'print ''    </TR>'';\n'];
     text = [text 'print ''    <TR>'';\n'];
     text = [text 'print ''    <TD align="center" width="800" bgcolor="#006699" height="23"><BR></TD>'';\n'];
     text = [text 'print ''    </TR>'';\n'];
     text = [text 'print ''  </TBODY>'';\n'];
     text = [text 'print ''<TABLE width="900" border="0" cellspacing="0" cellpadding="0" height="500">'';\n'];
     text = [text 'print ''  <TBODY>'';\n'];   
     text = [text 'print ''<TD rowspan="10" colspan="1" valign="top" align="center" width="10"></TD>'';\n'];
     text = [text 'print ''      <TD rowspan="4" colspan="2" valign="top" align="center" width="450" height="350">'';\n'];
     text = [text 'print ''      <TABLE width="300" ><TR><TD><BR></TD></TR></TABLE>'';\n'];
     text = [text 'print ''      <TABLE width="740" valign="top" cellspacing="1" cellpadding="4">'';\n'];
     text = [text 'print ''        <TBODY>'';\n'];
     text = [text 'print "<TR><TD bgcolor=\\"#DDDDDD\\" align=\\"center\\" width=\\"75\\" colspan=5><p class=blue>Gene name: $gene<br>No. of images in this clone: $nimg{$gene}</p></TD></TR>";\n'];
     text = [text 'print ''<TR><TD align="center" width="75" colspan=5><p class=blue></p></TD></TR>'';\n'];
     fprintf(fid, text);
     text = 'my $i = 1;\n';
     text = [text 'my $j = 1;\n'];
     text = [text 'while($i <= $quan) {\n'];
     text = [text '\tif ($j == 1) {\n'];
     text = [text '\t\tprint ''<TR>'';\n\t}\n'];
     text = [text '\tprint "<TD width=200 align=center><a href=\\"'...
             '$imgs[$i - 1]\\"><IMG src=\\"' ...
	     '$imgs[$i - 1]\\" width=180 height=144></a></TD>";\n'];
     text = [text '\t$i++;\n'];
     text = [text '\t$j++;\n'];
     text = [text '\tif ($j > 5) {\n'];
     text = [text '\t\t$j = 1;\n'];
     text = [text '\t\tprint ''</TR>'';\n\t}\n'];
     text = [text '};\n'];
     text = [text 'if ($j <= 5) {\n'];
     text = [text '\tprint ''</TR>'';\n}\n'];
     %text = [text 'print '''';\n'];
     fprintf(fid, text);

     text = 'print ''        </TBODY>'';\n';
     text = [text 'print ''      </TABLE>      </TD>    </TR>    </TR>'';\n'];
     text = [text 'print ''      <TR>        <TD colspan="2" align="center" width="450"><BR></TD>      </TR>    </TBODY>'';\n'];
     text = [text 'print ''  </TABLE>  </CENTER>'';\n'];
     text = [text 'print end_html;\n'];
     fprintf(fid, text);
     fclose(fid);
     unix(['chmod 755 ' cgifile]);
