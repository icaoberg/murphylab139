function [masks,mask_all] = ml_voronoi( dna,resolution,min_nucls_diam,max_nucls_diam,min_roundness1,min_roundness2)

%function [masks,mask_all] = ml_voronoi( dna,resolution,min_nucls_diameter,max_nucls_diam,min_roundness1)
% dna is the image matrix of the dna channel
% resolution: unit is micron per pixel
% min_nucls_diam: is the minimum diameter allowed for a nucleus, unit: micrometer
% max_nucls_diam: is the maxmum diameter allowed for a nucleus,unit micrometer
% min_roundness1: the threshold for pc, which describes the roundness (area of the object over area of convexhull around the object ), 50 was used by default
%% min_roundness2: the threshold for another roundness, which describes the relationship of perimeter and area. 70 was used by default
% Ivan, Aug4, 2005
% Modified by Yanhua Hu,using standard ml functions.

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

[size_x,size_y] = size(dna);

if ~exist('min_nucls_diam','var')
     min_nucls_diam=5; %micrometer
end

if ~exist('max_nucls_diam','var')
     max_nucls_diam=20; %micrometer
end 

if ~exist('min_roundness1','var')
     min_roundness1 = 50;
end

if ~exist('min_roundness2','var')
     min_roundness2 = 0.7;
end


masks={};
mask_all = [];
min_obj= (min_nucls_diam/resolution)^2;
max_obj= (max_nucls_diam/resolution)^2;

if length(unique(dna))==1
    return;
end

[dna,dnaBin] = ml_preprocess(double(dna),[],'ml','quantile');

imagelabeled = bwlabel(dnaBin) ;
obj_number = max(imagelabeled(:)) ;
if (obj_number==0)
   return
end
obj_sizes=[];
objx=[];
objy=[];
objCofsAll=[];
[img_moment00, img_moment10, img_moment01, obj_size] = ml_moments_1(int32(dna),int32(imagelabeled));

for (i=1:obj_number)
        obj_sizes(i) = obj_size(i+1);
        obj_m00 = double(img_moment00(i+1));
        obj_m10 = double(img_moment10(i+1));
        obj_m01 = double(img_moment01(i+1));
	objy(i)= obj_m10/obj_m00;
	objx(i)= obj_m01/obj_m00; 
end
objCofsAll_orig = [objx;objy];

% filter out the small noises here. The rest filtering is done after voronoi
objCofsAll=objCofsAll_orig(:,find(obj_sizes>=min_obj));
objx = objCofsAll(1,:);
objy = objCofsAll(2,:);

% object centers of the true nuclei 
objCofs = filter_nuclei( dna,dnaBin,obj_sizes,objCofsAll_orig,min_obj,max_obj,min_roundness1,min_roundness2);

%keyboard;

if isempty(objCofs)
  return;
end

number_of_nuclei = size( objCofs,2);

% Find the center of fluorescence of each of the objects
%[totalCof, objCofs] = ml_findCOFs( nuclei_obj, dna);

% find the sementation by voronoi for the all objects
%YH: add the four vertices as 4 seeds, so the code can handle case with less than 5 cells.
objx=[objx,1,1,size_x,size_x];
objy= [objy,1,size_y,1,size_y];
%[vx, vy] = voronoi( objx, objy );
failed = 0;
eval ('[vx, vy] = voronoi( objx, objy, {''QJ''} );', 'failed = 1');
if failed
    return;
end
newIm = zeros([size_x,size_y]);
%Connects all the vertices of the Voronoi diagram
for i = 1:size(vx,2)
    newIm = ml_draw_line2(vx(1,i), vy(1, i), vx(2, i), vy(2, i), newIm, 30);
end
%Connects the four corners of the image
newIm = ml_draw_line2(1,1,1,size_y,newIm, 30);
newIm = ml_draw_line2(1,size_y,size_x,size_y,newIm, 30);
newIm = ml_draw_line2(size_x,size_y,size_x,1, newIm, 30);
newIm = ml_draw_line2(size_x,1,1,1,newIm, 30);

% Creates binary imaged
segLines = im2bw(newIm, 1);
segLabel = bwlabel(imcomplement(segLines),4);

list = [];
dup = []; %duplicates, one mask contains more than one objcof
for i = 1:size(objCofs,2)
    r = segLabel(round(objCofs(1,i)), round(objCofs(2,i)));
    if ( find(list==r) )
        dup = [dup r];
    else
        list = [list r];
    end
end

mask_all=zeros([size_x,size_y]);
for i = 1 : number_of_nuclei
    mask = zeros([size_x,size_y]);
    % get nuclei region number;
    reg = segLabel(round(objCofs(1,i)), round(objCofs(2,i)));
    if (~sum(find(dup==reg)))
          mask(find(segLabel == reg)) = 1;
    end
    masks{i} = uint8(mask);
    mask_all=mask_all+mask;
end




%%%%%%%% helper function
 
function objCofs = filter_nuclei( dna,dnaBin,obj_sizes,objCofsAll,min_obj,max_obj,min_roundness1, min_roundness2);

rng_fac =2;
%dnaBin=im2bw(dna);

bw_Perimeter = bwperim(dnaBin);
label_Obj = bwlabel(dnaBin);
label_perimeter= label_Obj;
label_perimeter(find(bw_Perimeter==0))=0;

objLabel = unique(label_Obj);
num_objs = (length(objLabel)-1);
objPerims = [];
edgetouch = zeros(1,num_objs);
imgsize=size(dna);

% Figure out which ones touch the edge
% the roundness1,
% And perimeter of each objects

for i = 1 : num_objs
    [xs,ys]=find(label_Obj==i);
    isEdge=sum(xs==1)+sum(xs==imgsize(1))+sum(ys==1)+sum(ys==imgsize(2));
    if isEdge
       edgetouch(i) = 1;
    end
 
    % if it is a line, convhull cannot be calculated
    if length(unique(xs))<4|length(unique(ys))<4
        roundness1(i)=0;
        roundness2(i)=0;
        continue
    end
    [k,ar] = convhull(xs',ys');
    if ar > 0 
         roundness1(i) = round(obj_sizes(i)*100/ar);
    else
        roundness1(i) = 0;
    end
  
    objPerim=length(find(label_perimeter==i));
    roundness2(i) = 4*pi*obj_sizes(i)./(objPerim.^2)
   
end




obj_sizes_Qual = obj_sizes(find(obj_sizes>=min_obj));
med = mean(obj_sizes_Qual);
std_size = std(obj_sizes_Qual);
rng = round(rng_fac*std_size);
min_limit = max(med-rng,min_obj);
%min_limit = med-rng;
max_limit = min(med+rng,max_obj);

%keyboard;
% Keep only nucleus candidates which are within range 

goodnuclei = [ find(edgetouch==0 & obj_sizes>=min_limit & obj_sizes<= max_limit  & roundness1 >= min_roundness1 & roundness2>= min_roundness2 )];

%goodnuclei = [ find(obj_sizes>=(med-rng) & obj_sizes<=(med+rng) & obj_sizes<=max_obj & roundness1 >= min_roundness1 )];


objCofs = objCofsAll(:,goodnuclei);
