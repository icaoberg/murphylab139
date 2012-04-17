function cocmat = ml_updtCOCMAT(img1,img2,cocmat,graylevel) 

%ML_UPDTCOCMAT update the existing temporal co-occurence matrix 'cocmat' using img1 and img2
%   COCMAT = ML_UPDTCOCMAT(IMG1,IMG2,COCMAT,GRAYLEVEL)

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

%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%  
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%  
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%  
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu

%   Yanhua Hu, Mar,05

s=length(img1(:));
if isempty(cocmat)
     cocmat=zeros(graylevel);
end


for i=1:s
  a=img1(i);
  b=img2(i);
  if a~=0&b~=0
     cocmat(a,b)= cocmat(a,b)+1; 
     cocmat(b,a)=cocmat(b,a)+1;
  end
end
   

