%This script loads UCSF labeling information from allOrfData.txt into a variable yeaststruct (this may takes about two minutes)

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

filename = ['data' filesep 'allOrfData.txt'];
nmax = 32;
outputstr = '[';
ms = '';
for n = 1:nmax
%    eval(['M' num2str(n) ' = magic(n)']);
    if  n == nmax
        sshort = sprintf('s{%d}]',n);          
    else 
        sshort = sprintf('s{%d},',n);    
    end
    outputstr = strcat(outputstr,sshort);
    ms = strcat(ms,'%s');
end

textreadcmd = sprintf('%s = textread(''%s'',''%s%%*[^\\n]'',''delimiter'',''\\t'');', outputstr, filename, ms);
eval(textreadcmd);

for n = 1:nmax
    cmd = sprintf('s{%d}{1}(find(s{%d}{1} == '' '')) = ''_'';', n, n);  eval(cmd);
    cmd = sprintf('s{%d}{1}(find(s{%d}{1} == ''?'')) = '''';', n, n);  eval(cmd);
end

lmax = length(s{1});
% the 9th cell is the location information
for l = 2:lmax
    if ~isempty(find(s{9}{l} == ','))
        s{9}{l} = s{9}{l}(1:find(s{9}{l} == ',')-1);
    end
end

for l = 2:lmax
    for n  = 1:nmax
        sshort = eval(sprintf('s{%d}{1},',n)); 
        cmd = sprintf('yeaststruct.%s{%d} = s{%d}{%d};',sshort,l-1,n,l);
        eval([cmd]);
    end
end

