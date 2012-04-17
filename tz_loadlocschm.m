function loc = tz_loadlocschm(filepath)
%TZ_LOADLOCSCHM Load localization schema from a file.
%   LOC = TZ_LOADLOCSCHM(FILEPATH) returns an array of structures. Each
%   elements have two fields:
%       'label' - number label
%       'name' - location name
%   
%   See also

%   04-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end

str = textread(filepath,'%s','delimiter','\n');

for k=1:length(str)
    [loc(k).label,name] = strtok(str{k});
    loc(k).name = strtrim(name);
end
