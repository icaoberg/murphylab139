function [lines, masks] = ml_segment_cells_voronoi_super(im, dnaIm)

%function [lines, masks] = ml_segment_cells_voronoi_super(im, dnaIm)
% im is the protein image
% dnaIm is the image of the dna channel
% Lines is the image with borders to seperate each cell\
% masks is a cell array containing masks

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
     dna = double(dnaIm);
     dna = ml_imgbgsub(dna, 'common');
    
     
     % Threshold
     threshDna = 0;
     eval('threshDna = ml_choosethresh(dna);', ...
	  'disp(''could not get threshold''); threshDna = 10');
     dnaBin = im2bw(uint8(dna), threshDna/255); %65535
      
%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% OBJECTFINDING %%%%%%%%%%%%%%%%%%%
     
     findholes = 0;
     dnaobj = ml_3dfindobj( dnaBin, findholes);
     si = uint8(split_im(dnaBin));
     
     save('si', 'si');
    % disp('End blah');
    % keyboard;
     % Find those objects that are large enough to be nuclei and are
     % not touching the edge of the image
     %nuclei_obj = filter_nuclei( dnaobj, size(dnaBin));
     nuclei_obj = filter_nuclei2(dnaobj, size(si), si);
     
%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% SEEDIM %%%%%%%%%%%%%%%%%%%%%%%
  
     seedimg = uint8(zeros(size(si)));
     number_of_nuclei = length(nuclei_obj);
     for i = 1 : number_of_nuclei
	 v = double(nuclei_obj{i}.voxels);
	 seedimg(sub2ind(size(si),v(1,:),v(2,:))) = i;
     end
     % Add an additional seed region containing the edges of the image
     % Thiims is to help get rid of partial cells
     seedimg(:,[1 end]) = number_of_nuclei + 1;
     seedimg([1 end],:) = number_of_nuclei + 1;
     %imshow(seedimg, []);

%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%


     
%%%%%%%%%%%%%%%% LOCATION FINDER %%%%%%%%%%%%%%%%%%
     addpath /imaging2/backup/imaging/usr/matlab6/toolbox/images/images/
     [totalCof, objCofs] = ml_findCOFs( nuclei_obj, dna);
     siObjs = ml_3dfindobj(si, 0);
     [totalCofAll, objCofsAll] = ml_findCOFs(siObjs, dna);
     %keyboard;
%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
    
     
     objy = objCofsAll(1,:);
     objx = objCofsAll(2,:);
    
     %bw1 = bwmorph(dnaBin, 'thin', Inf);
     %[bi,bj] = ind2sub(size(bw1), find(bw1));
     %bi = abs((bi-1)-size(bw1,2));
     %bj = abs((bj-1)-size(bw1,1));
     %objx = bi'; 
     %objy = bj';
     
     s = size(objy);
     objx = [objx,0,0,size(dna,2),size(dna,2)];
     objy = [objy,0,size(dna,1), 0 ,size(dna,1)];
%%%%%%%%%%%%%%%%%%%% DELAUNAY %%%%%%%%%%%%%%%%%%%%%%%%
     fprintf('Delaunay TRI\n');
     TRI = delaunay(objx, objy);

     [vx, vy] = voronoi(objx, objy, TRI);
     newIm = zeros(size(dna));
     fprintf('Drawing %.0f lines\n', size(vx,2));
     for i = 1:size(vx,2)
	 newIm = ml_draw_line(vy(1,i), vx(1, i), vy(2, i), vx(2, i), ...
			      newIm, 30);
	 fprintf('.');
     end
     fprintf('\n');
     newIm = ml_draw_line(1,1,1,size(newIm,2),newIm, 30);
     newIm = ml_draw_line(1,size(newIm,2),size(newIm,1), ...
                           size(newIm,2),newIm, 30);
     newIm = ml_draw_line(size(newIm,1),size(newIm,2), ...
                           size(newIm,1),1, newIm, 30);
     newIm = ml_draw_line(size(newIm,1),1,1,1,newIm, 30);
     
   
     segLines = im2bw(newIm, 1);
     lines = segLines;
%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %keyboard;
     

%%%%%%%%%%%%%%%%% MASK %%%%%%%%%%%%%%%%%%%%%%%%%%%

    segLabel = bwlabel(imcomplement(segLines),4);
    
    list = [];
    exc = [];
    for i = 1:size(objCofs,2)
	r = segLabel(round(objCofs(1,i)), round(objCofs(2,i)));
	if (~isempty(list)) 
	    if (find(list==r))
		exc = [exc r];
	    else
		list = [list r];
	    end
	else
	    list = [list r];
	end
	
    end
    
    
    %for i = 1 : number_of_nuclei
    %	mask = zeros(size(dna));
    %	% get nuclei region number;
    %	reg = segLabel(round(objCofs(1,i)), round(objCofs(2,i)));
    %	if (~isempty(exc))
    %	    if (find(exc~=reg))
    %		mask(find(segLabel == reg)) = 1;
    %	    end
    %	end
    %	masks{i} = mask;
    %    end
    masks = {};
    for i = 1:length(list)
	if (isempty(exc) | isempty(find(exc==i)))
	    mask = zeros(size(dna));
	    mask(find(segLabel == list(i))) = 1;
	    masks{length(masks)+1} = mask;
	    %ml_superimpose(mask, dnaBin); pause;
	end
    end
    
	
%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%




function nucobj = filter_nuclei( obj, imgsize)
% Find the objects in the OBJ list that are likely to be nuclei
% and are not touching the edge
 
% Get a list of all object sizes
for i = 1 : length(obj)
    objsizes(i) = obj{i}.size;
end

% Figure out which ones touch the edge
edgetouch = ones(1,length(obj));
im = uint16(zeros(imgsize));
im2 = double(zeros(imgsize));
%im2 = im2*65535;
cm = [0 0 0; 1 1 1; 1 0 0; 0 1 1];

for oi = 1 : length(obj)
    o = obj{oi};
    if( find(o.voxels([1 2],:)==1))
    else
	if( find(o.voxels(1,:)==imgsize(1)))
	else
	    if( find(o.voxels(2,:)==imgsize(2)))
	    else
		edgetouch(oi) = 0;
	    end
	end
    end
    k = convhull(double(o.voxels(2,:)),double(o.voxels(1,:)));
    ar = polyarea(double(o.voxels(2,k)),double(o.voxels(1,k)));
    
    %im2 = (im2~=0);
    %v = double(o.voxels);
    %im2(sub2ind(imgsize, v(1,:), v(2,:))) = 1;
    %end
    
    if ar > 0 
	pc(oi) = round((o.size)*100/ar);
    else
	pc(i) = 100;
    end
    
end

%%%%%%%%%%% OBJECT SIZE FILTER %%%%%%%%%%%%
med = median(objsizes);
rng = round(med/3);
%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%
% Keep only nucleus candidates which are within range and do not
% touch the edge
goodnuclei = [ find(edgetouch==0 & objsizes>(med-rng) & objsizes<(med+rng) ...
		    & pc >= 99 )];
nucobj = obj(goodnuclei);





function nucobj = filter_nuclei2(obj, imgsize, im)

% Get a list of all object sizes
for i = 1:length(obj)
    objsizes(i) = obj{i}.size;
end

%Figure out which ones touch the edge
edgetouch = ones(1,length(obj));
multiDna = ones(1,length(obj));
wrongSize = ones(1,length(obj));
cm = [0 0 0; 1 0 0; 1 1 1; 0 1 0; 0 0 1];


im2 = uint8(zeros(imgsize));


for oi = 1:length(obj)
    o = obj{oi};
    
    %EDGETOUCH FILTER
    for i = 1:o.size
	im2(o.voxels(1,i), o.voxels(2,i)) = 1;
    end
    if (find(o.voxels([1 2],:)==1))
    else
	if (find(o.voxels(1,:)==imgsize(1)))
	else
	    if(find(o.voxels(2,:)==imgsize(2)))
	    else
		edgetouch(oi) = 0;
	    end
	end
    end
   
    
    %MULTIDNA FILTER
    k = convhull(double(o.voxels(2,:)), double(o.voxels(1,:)));
    tempIm = zeros(size(im2));
    if size(k,2) > 1
	for i = 1:length(k)-1
	    st = double([o.voxels(1,k(i)), o.voxels(2,k(i))]);
	    en = double([o.voxels(1,k(i+1)), o.voxels(2,k(i+1))]);
	    tempIm = ml_draw_line(st(1), st(2), en(1), en(2), ...
				   tempIm, 3);
	end
	addArea = length(find(tempIm==3));
	ar = floor(polyarea(double(o.voxels(2,k)),double(o.voxels(1,k))));
	ar = ar + addArea;
    else
	ar = 1;
	addArea = 1;
	tempIm(double(o.voxels(1,1)),double(o.voxels(2,1))) = 3;
    end
	%fprintf('Convhull Area: %3.0f\n', ar);
	%fprintf('Convhull without addition: %3.0f\n', ar-addArea);
	%fprintf('Actual Area: %3.0f\n', double(o.size));
	%fprintf('Percent diff: %3.2f\n', double(o.size)*100/ar);
    singleObj = zeros(size(im2));
    singleObj(find(im2 == 1)) = 1;
    singleObj(find(tempIm==3)) = 1;
    h = ml_3dfindobj(uint8(singleObj), 1);
    
    if (length(h) ~= 1) error('Too many objects'); end
    tempSize = 0;
    if (h{1}.n_holes > 0)
	for j = 1:h{1}.n_holes
	    if (h{1}.holes{j}.size > 1)
		tempSize = tempSize + h{1}.holes{j}.size;
	    end
	end
    end
    
    if (round(tempSize*100/ar) < 1)
	multiDna(oi) = 0;
    end
    
	%fprintf('Hole size total: %3.0f\n', tempSize);
	fprintf('Percent Holiness: %3.0f\n', round(tempSize*100/ar));
	%fprintf('Percent of Object in Con Hull Line: %3.2f\n\n', ...
	
	%Add code here to recursively reduce objects until
	% there are no more holes in a convex hull binded box
	%A single object should create more objects corresponding
	% to single nuclei
	%However, you may just want to discard touching/overlapping
	% DNA. In which case...add nothing
	
	%fprintf('Number of overlapping pixels: %3.0f\n', ...
	%	length(intersect(find(im2==1), find(tempIm==3))));
	%fprintf('Number of pixels in hull line: %3.0f\n', ...
	%	addArea);
	%fprintf('Percent of hull pixels in image: %3.0f\n', ...
	%	length(intersect(find(im2==1), find(tempIm==3)))* ...
	%	100/addArea);
	%singleObj = zeros(size(im2));
	%singleObj(find(im2 == 1)) = 1;
	%while (round(tempSize*100/ar) > 0)      
	%    singleObj = mmdist(singleObj);
	%    singleObj = uint16(double(singleObj) -1);
	%    singleObj(find(singleObj<0)) = 0;
	%   imshow(singleObj, []); pause; 
	%    singleObj(find(tempIm==3)) = 1;
	%    h = ml_3dfindobj(uint8(singleObj), 1);
	%    %if oi == 12 keyboard; end
	%    if (length(h) ~= 1) error('Too many objects'); end
	%    tempSize = 0;
	%    if (h{1}.n_holes > 0)
	%	for j = 1:h{1}.n_holes
	%	    if (h{1}.holes{j}.size > 1)
	%		tempSize = tempSize + h{1}.holes{j}.size;
	%	    end
	%	end
	%    end
	%    fprintf('Hole size total: %3.0f\n', tempSize);
	%    fprintf('Percent Holiness: %3.0f\n', tempSize*100/ar);
	%   
	%    
	%end
    %SIZE FILTER
    if (objsizes(oi) > round(median(objsizes) - median(objsizes)/3) ...
	& objsizes(oi < round(median(objsizes) + median(objsizes)/3)))
	wrongSize(oi) = 0;
    end
    
    im2 = uint16(double(im2)+double(tempIm));
    % To look at each object checked, uncomment the following line
    % db = 1; 
    
    if (exist('db'))
	inc = ((edgetouch(oi) + multiDna(oi) + wrongSize(oi)) == 0);
	if (inc)
	    tempSt = 'included';
	else
	    tempSt = 'excluded due to:';
	end
	
	fprintf('  Object %.0f %s %s%s%s\n', oi, tempSt, ...
		char('|touching edge|' * edgetouch(oi)), ...
		char(['|touching DNA (hole pct is ' ...
		     num2str(round(tempSize*100/ar)) ')|'] ...
		     * multiDna(oi)), ...
		char('|wrong size|' * wrongSize(oi)));
	imshow(im2, cm);  pause; close;
	
	
    end
    im2(find(im2 == 1)) = 2;
end

goodnuclei = [find(edgetouch==0 & multiDna==0 & wrongSize==0)];
nucobj = obj(goodnuclei);

function split_im = split_im(im)
%Find objects in image
initobj = ml_3dfindobj( im, 0);
splitObjs = {};
split_im = zeros(size(im));
objs = {};
objsToDo = initobj;
for i = 1:length(initobj)
    if (~isstruct(initobj{i}))
	disp('Missing init obj');
	keyboard;
    end
end
%For each object
while(~isempty(objsToDo))
    fprintf('To do: %.0f\n', length(objsToDo));
    o = objsToDo{1};
    if (~isstruct(o))
	keyboard;
    end
    
    oIm = zeros(size(im));
  
    for i = 1:o.size
	oIm(o.voxels(1,i), o.voxels(2,i)) = 1;
    end
       
    pct = convHoles(o, size(oIm));
    
    if (round(pct) > 1)
	%	%If 1% of the image are holes, it's a multiDNA object,
	%	% so split the object
        %	disp('Split time');
		objs = split_obj(o, size(oIm));
    else
	splitObjs{length(splitObjs)+1} = o;
    end
    
    if length(objsToDo) == 1
	objsToDo = {};
    else
	objsToDo = {objsToDo{2:end}}';
    end
    
    if length(objs)<=1
	continue;
    end
    
    preLength = length(objsToDo);
    for i = 1:length(objs)
	%if length(splitObjs) == 
	%    disp('Error alert');
	%    keyboard;
	%end
	%for k = 1:length(objsToDo)
	%    if (~isstruct(objsToDo{k}))
	%	disp('Error previously assigning to objsToDo');
	%	keyboard;
	%    end
	%end
	if (~isstruct(objs{i}))
	    disp('huh?');
	    keyboard;
	end
	
	objsToDo{preLength+i} = objs{i};
	for k = 1:length(objsToDo)
	    if (~isstruct(objsToDo{k}))
	    disp('Error assigning to objsToDo');
	    keyboard;
	    end
	end
	
    end
end

for i = 1:length(splitObjs)
    for j = 1:splitObjs{i}.size
	o = splitObjs{i};
	split_im(o.voxels(1,j), o.voxels(2,j)) = 1;
    end
end


function objs = split_obj(o, im_size)
oIm = uint8(zeros(im_size));
for i = 1:o.size
    oIm(o.voxels(1,i), o.voxels(2,i)) = 1;
end

de = mmdist(im2bw(oIm, 1/255), mmsebox, 'EUCLIDEAN');

numOfObjs = 1;
thresh = 0;
while( numOfObjs <= 1)
   ind = [];
   thresh = thresh + 1;
   ind = find(de>thresh);
   if (length(ind)==0)
       objs = {o};
       break;
   end
   de2 = zeros(size(de));
   de2(ind) = de(ind);
   de = de2;
   
   newObjs = ml_3dfindobj(uint8(im2bw(de, 1/255)), 0);
   numOfObjs = length(newObjs);
end

de = im2bw(de, 1/255);
oIm = roifilt2(0, oIm, ~de);
objs = ml_3dfindobj(oIm, 0);
for i = 1:length(objs)
    if (~isstruct(objs{i}))
	keyboard;
    end
end


function pct = convHoles(o, im_size)
cm = [0 0 0; 1 0 0; 1 1 1; 0 1 0; 0 0 1];
% Object: 1
% Conv line: 3

% Put the object in an image
oIm = uint8(zeros(im_size));
for i = 1:o.size
    oIm(o.voxels(1,i), o.voxels(2,i)) = 1;
end

% Calculate convex hull of object
k = convhull(double(o.voxels(2,:)), double(o.voxels(1,:)));

% Put convex hull in an image
convIm = uint8(zeros(im_size));
if size(k,2) > 1
    for i = 1:length(k)-1
	st = double([o.voxels(1,k(i)), o.voxels(2,k(i))]);
	en = double([o.voxels(1,k(i+1)), o.voxels(2,k(i+1))]);
	convIm = ml_draw_line(st(1), st(2), en(1), en(2), ...
			      convIm, 3);
    end
    addArea = length(find(convIm==3));
    ar = floor(polyarea(double(o.voxels(2,k)),double(o.voxels(1,k))));
    ar = ar + addArea;
else
    ar = 1;
    addArea = 1;
    convIm(double(o.voxels(1,1)),double(o.voxels(2,1))) = 3;
end

% Combine object image and convhull image
singleObj = uint8(zeros(im_size));
singleObj(find(oIm == 1)) = 1;
singleObj(find(convIm==3)) = 1;

%imshow(uint8(double(oIm)+double(convIm)), cm); keyboard;
% Find holes. H here is the object containing the holes.
%  It must be 1.
h = ml_3dfindobj(uint8(singleObj), 1);
if (length(h) ~= 1) 
    warning('Not 1 object');
    pct = NaN;
    return;
end


holesSize = 0;
if (h{1}.n_holes > 0)
    for j = 1:h{1}.n_holes
	if (h{1}.holes{j}.size > 1)
	    holesSize = holesSize + h{1}.holes{j}.size;
	end
    end
end

pct = holesSize*100/ar;
