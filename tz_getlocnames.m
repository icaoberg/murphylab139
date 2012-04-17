function locnames = tz_getlocnames(locschms,loclabels)
%TZ_GETLOCNAMES Get [location name]s from location schemes.
%   LOCNAMES = TZ_GETLOCNAMES(LOCS,LOCLABELS) returns a cell array of
%   [location name]s of the cell array of [locaion label]s LOCLABELS.
%   LOCSCHMS provides a map between them.
%   
%   See also

%   06-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

locschms = squeeze(struct2cell(locschms));

locindices = strmatch_multi(loclabels,locschms(1,:));

alllocnames = {locschms{2,:} 'unknown'};
locindices(find(locindices==0)) = length(alllocnames);
locnames = alllocnames(locindices);