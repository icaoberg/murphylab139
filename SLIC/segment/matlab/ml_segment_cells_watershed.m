function [lines, segment_masks] = ml_segment_cells_watershed( cellimg, dnaimg)

% SEGMENT_MASKS = ML_SEGMENT_CELLS_WATERSHED( CELLIMG, DNAIMG)
%
% Segment cells in the CELLIMG using the watershed algorith
% and using the nuclei from the DNA image to seed the watershed.

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


%%%%%%%%%%%% PROCESSING %%%%%%%%%%%%%%%%%%%%%%
     cellimg = double(cellimg);
     dna = double(dnaimg);
     dna = ml_imgbgsub(dna, 'common');
    
     
     % Threshold
     threshDna = 0;
     eval('threshDna = ml_choosethresh(dna);', ...
	  'disp(''could not get threshold''); threshDna = 10');
     dnaBin = im2bw(uint8(dna), threshDna/255);
     
%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% OBJECTFINDING %%%%%%%%%%%%%%%%%%%
     
     findholes = 0;
     dnaobj = ml_3dfindobj( dnaBin, findholes);
     % Find those objects that are large enough to be nuclei and are
     % not touching the edge of the image
     nuclei_obj = filter_nuclei( dnaobj, size(dnaBin));
     

%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% SEEDIM %%%%%%%%%%%%%%%%%%%%%%%
  
     seedimg = uint8(zeros(size(dnaBin)));
     number_of_nuclei = length( nuclei_obj);
     for i = 1 : number_of_nuclei
	 v = double(nuclei_obj{i}.voxels);
	 seedimg(sub2ind(size(dnaBin),v(1,:),v(2,:))) = i;
     end
     % Add an additional seed region containing the edges of the image
     % This is to help get rid of partial cells
     seedimg(:,[1 end]) = number_of_nuclei + 1;
     seedimg([1 end],:) = number_of_nuclei + 1;

% Do the watershed
  inv_cellimg = uint16(max(cellimg(:))-cellimg);
  regions = 0;
  success = 0;
  while( success == 0)
      success = 1;
      eval('regions = mmcwatershed( inv_cellimg, seedimg);', ...
	   'disp(''Trying mmcwatershed again''); success=0;');
  end
  
  %regions = mmcwatershed( inv_cellimg, seedimg);

  lines = regions;
% Label the regions
  regions = 1 - double(regions);
  labels = bwlabel( regions, 4);
  edge_label = labels(1,1);
  labels(find(labels==edge_label))=0;
  labels(find(labels))=1;
  segment_masks = bwlabel( labels);
  
  
  
function nucobj = filter_nuclei( obj, imgsize)
% Find the objects in the OBJ list that are likely to be nuclei
% and are not touching the edge

% Get a list of all object sizes
for i = 1 : length(obj)
    objsizes(i) = obj{i}.size;
end

% Figure out which ones touch the edge
edgetouch = ones(1,length(obj));
for i = 1 : length(obj)
    o = obj{i};
    if( find(o.voxels([1 2],:)==1))
    else
	if( find(o.voxels(1,:)==imgsize(1)))
	else
	    if( find(o.voxels(2,:)==imgsize(2)))
	    else
		edgetouch(i) = 0;
	    end
	end
    end
end

% Keep only nucleus candidates which are big enough and do not
% touch the edge
M = max(objsizes);
goodnuclei = [ find(edgetouch==0 & objsizes>M/3)];
nucobj = obj(goodnuclei);
