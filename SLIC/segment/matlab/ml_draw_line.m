function I = ml_draw_line(r1, c1, r2, c2, I, lineVal)
%ML_DRAW_LINE Draw a line in an image.
%   I = ML_DRAW_LINE(R1,C1,R2,C2,I,LINEVAL) returns an image that is the
%   superimpose of the input image I and a line from [R1,C1] to [R2,C2].
%   LINEVAL is the intensity of the line.
%   
%   See also

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

insize = size(I);
if(r1==r2 & c1==c2)
    return;
end


if (abs(r2-r1) > abs(c2-c1))
    for i = min(r1,r2):max(r1,r2)
	%if (i < 1 | i > size(I,1)) continue; end;
	if (i < 1) 
	    if (max(r1,r2) < 1) return; else i = 1; end;
	else 
	    if (i > size(I,1)) return; end
	end
	c = interp1([r1 r2], [c1 c2], round(i), 'spline');
	if (c < 1 | c > size(I,2)) continue; end
	I(round(i), round(c)) = lineVal;
    end
else
    for i = min(c1,c2):max(c1,c2)
	%if (i < 1 | i > size(I,2)) continue; end;
	if (i < 1)
	    if (max(c1,c2) < 1) return; else i = 1; end;
	else 
	    if (i > size(I,2)) return; end
	end
	r = interp1([c1 c2], [r1 r2], round(i), 'spline');
	if (r < 1 | r > size(I,1)) continue; end
	I(round(r), round(i)) = lineVal;
    end
end


