function [t,loclabel,mainlabels] = tz_islocunique(loclabel)
%TZ_ISLOCUNIQUE Extract unique [location label]s
%   T = TZ_ISLOCUNIQUE(LOCLABEL) returns the 
%   
%   See also

%   04-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

if isempty(loclabel)
    t = 0;
    return;
end

mainlabels = strtok(loclabel,'.');
loclabel = unique(mainlabels);
t = length(loclabel)==1;

