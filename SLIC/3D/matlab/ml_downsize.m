function new = ml_downsize(old, ratio, dsmeth)
% FUNCTION NEW = ML_DOWNSIZE(OLD, RATIO, DSMETH)
% Downsize a 3D image using the downsize ratio
% Important note: this version just works for int ratio value.  Check
%     ml_3dimresize for real ratios
% old: the original image
% ratio: a length 3 vector, ratio(1) is the downsample ratio on x, ratio(2)
%        for y and ratio(3) for z.  All 3 dimensions will be halfed if a
%        [2 2 2] is used as ratio
% dsmeth: sum ('summation',blank) or average ('average'), jnewberg 11/24/05
% 
% Xiang Chen
% Aug 15, 2002

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

if ~exist( 'dsmeth','var')
    dsmeth = 'summation';
end


SIZE = size(old);
ratio = round(ratio);
ratio(find(ratio == 0)) = 1;
if ~(ratio - [1 1 1])
    new = old;
    return;
else
    new = zeros(floor(SIZE(1) / ratio(1)), ...
		floor(SIZE(2) / ratio(2)), ...
		floor(SIZE(3) / ratio(3)));
end

for m = 1 : floor(SIZE(3)/ratio(3))
    for n = 1 : floor(SIZE(2) / ratio(2))
        for o = 1 : floor(SIZE(1) / ratio(1))
            tmp = old(1 + (o - 1) * ratio(1) : o * ratio(1), ...
		      1 + (n - 1) * ratio(2) : n * ratio(2), ...
		      1 + (m - 1) * ratio(3) : m * ratio(3));
            if ~strcmp(dsmeth,'average')
                new(o, n, m) = sum(tmp(:));
            else
                % added Xiang's average method, jnewberg, 11/24/05
                new(o, n, m) = uint8(round(sum(tmp(:)) / (ratio(1) * ratio(2) * ratio(3))));
            end
        end
    end
end
