function threshold = ml_rcthreshold(image)  
% ML_RCTHRESHOLD - Find a threshold for IMAGE 
% ML_RCTHRESHOLD(IMAGE) - input is output from nihscale(image)
% output is threshold value using same algorithm as nih image
% Ridler and Calvard (IEEE tans. on Systems, man, cybernetics
% SMC-8 no 8 aug 1978)
% Trussel (same journal) SMC-9 no 5 may 1979
%
% Basically, the algorithm moves the proposed threshold from the 
% minimum non-zero bin of the histogram to the maximum non-zero bin.
% While doing so, it calculates the 'center of mass' separately
% for that part of the histogram below the proposed threshold and
% that part above.  While the proposed threshold is less than 
% the mean of these two centers of mass, the algorithm continues.
% When this condition no longer holds, the proposed threshold is
% rounded to give a final result.
%
% If <= one intensity value is nonzero output threshold will be 1/2
% and error message will be printed to screen
%
% 20 Aug 98 - M.V. Boland.  Modified from version written by W. Dirks.
%
% Sometime in Year 2000, MV and RFM noticed that the algorithm
% had not been implemented properly. Corrected by MV and RFM.

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

% $Id: mb_nihthreshold.m,v 1.2 1999/02/17 14:19:57 boland Exp $

%
% Generate a 256 bin histogram
%
histo = imhist(image(:),256);

hsum = sum(histo);


histo(1) = 0;
histo(256) = 0;
error = 0; 

MinIndex = 1;
MaxIndex = 256;

while ( (histo(MinIndex) == 0) & (MinIndex < 256) )
	MinIndex = MinIndex + 1;
end;

while ((histo(MaxIndex) == 0) & (MaxIndex > 1))
	MaxIndex = MaxIndex - 1;
end;

if (MinIndex >= MaxIndex)
	Result = 128;
	disp('There is <= 1 nonzero pixel intensity');
	error = 1;
end;

MovingIndex = MinIndex;
if (error == 0)

	Prev = MinIndex-1;
	while (Prev < MovingIndex)
		if (MovingIndex > (MaxIndex-1))
			break;
		end;

		sum1 = sum([MinIndex:MovingIndex]' .* ...
                            histo(MinIndex:MovingIndex)) ;
		sum2 = sum(histo(MinIndex:MovingIndex)) ;
		sum3 = sum([MovingIndex+1:MaxIndex]' .* ...
                            histo(MovingIndex+1:MaxIndex)) ;
		sum4 = sum(histo(MovingIndex+1:MaxIndex)) ;

		Prev = MovingIndex;
		MovingIndex = floor((sum1/sum2 + sum3/sum4)/2) ;
                         
loweravg=sum1/sum2;
higheravg=sum3/sum4;

	end;

end;

threshold = MovingIndex - 1;
