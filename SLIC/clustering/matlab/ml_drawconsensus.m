function label_o = ml_drawconsensus(linksin, namein, fontsize, method, linewidth)
% FUNCTION ML_DRAWCONSENSUS(LINKSIN, NAMEIN, FONTSIZE, METHOD, LINEWIDTH)
% This function draws the consensus tree based on the result from 
% ml_readconsensus.
% linksin: links object from ml_readconsensus
% namein: names of the clones
% fontsize: font size
% method: 0 for Euclidean (default) and 1 for Mahalanobis
% linewidth: line width
% label_o: sequence of the labels on the graph, from left to right
% Xiang Chen, Jan 04, 2005
% XC, Aug 15, 2005, add output.

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

if (~exist('fontsize', 'var') | isempty(fontsize))
     fontsize = 10;
end

if (~exist('method'))
     method = 0;
end

if (~exist('linewidth', 'var') | isempty(linewidth))
     linewidth = 1;
end

links = linksin;
name = namein;
maxn = 0;
thresh = 35;
for m = 1 : length(name)
     if (length(name{m}) > maxn)
         maxn = length(name{m});
	 if (maxn > thresh)
             maxn = thresh;
	 end
     end
end


for m = 1 : length(name)
     if (length(name{m}) < maxn)
         name{m}(length(name{m}) + 1 : maxn) = ' ';
     end
end

             
ord = [];
for m = 1 : length(links)
     for n = 1 : length(links(m).Clusters)
         if (links(m).Clusters(n) <= length(name))
             ord = [ord links(m).Clusters(n)];
	     links(m).Clusters(n) = length(ord);
	 end
     end
end

maxd = links(end).Max;
startpos = maxd - xc_length(links);
label = name(ord');
label_o = namein(ord');
labelxoff = 0;
labelyoff = -0.1;
textrot = 270;

figure
set(gcf,'Position', [50, 100, 600, 300]);

X = 1:length(ord);
Y = startpos; 
A = [];
B = [];
for m = 1 : length(links)
    cluster = links(m).Clusters;
    d = links(m).Distance;
    d1 = d(1) + Y(cluster(1));
    for n = 1 : length(cluster)
        A = [A, [X(cluster(n)), X(cluster(n))]'];
	B = [B [Y(cluster(n)), d1]'];
        %plot([X(cluster(n)), X(cluster(n))], [Y(cluster(n)), d], 'k');
    end
    A = [A, [min(X(cluster)), max(X(cluster))]'];
    B = [B, [d1, d1]'];
    %plot([X(cluster(1)), X(cluster(n))], [d, d], 'k');
    X = [X mean(X(cluster))];
    Y = [Y d1];
end
axes('position', [0.05, 0.3, 0.92, 0.7]);

h = plot(A, B, 'k-');
for m = 1 : length(h)
    set(h(m),'LineWidth',linewidth);
end
set(gca,'FontSize',fontsize,'LineWidth',linewidth);
axis([0 length(ord)+2 0 maxd*1.05])
set(gca,'XTickLabel',[],'XTick',[],'box','off','xcolor','w');
h = text((1:length(ord))+labelxoff,zeros(length(ord),1)+labelyoff,label,'rotation',textrot, 'Fontsize', fontsize, 'Fontname', 'Courier');
if (method == 1)
     ylabel('Z-scored Mahalanobis Distance');
elseif (method == 2)
     ylabel('Z-scored City Block Distance');
elseif (method == 3)
     ylabel('Z-scored Cosine Angle Distance');

else
     ylabel('Z-scored Euclidean Distance');
end


function l = xc_length(links)
     ll = 1 - length(links);
     for m = 1 : length(links)
         ll = ll + length(links(m).Clusters);
     end
     tmp = zeros(1, length(links) + ll);
     for m = 0 : length(links) - 1
         x = length(links) - m;
         tmp(links(x).Clusters) = links(x).Distance + tmp(ll + x);
     end

     l = tmp(1:ll);
