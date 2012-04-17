function machine = tz_machineinfo(varargin)
%TZ_MACHINEINFO Get information of the machine
%   TZ_MACHINEINFO(INFO1,...,INFON) returns machine information, which is 
%   stored in a struecture.
%
%   Options for information:
%       'name': machine name
%       'computer': computer type (for details, type 'help computer')

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

%   24-May-2005 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

for i=1:length(varargin)
    switch varargin{i}
        case 'name'
            if exist('machine.mat','file')
                load machine.mat
            else
                [s,hostname]=unix('hostname');
                if uint8(hostname(end))==10
                    hostname=hostname(1:end-1);
                end
                dotpos=find(hostname=='.');
                if isempty(dotpos)
                    machine.name=hostname;
                else
                    machine.name=hostname(1:dotpos-1);
                end
            end
        case 'computer'
            machine.computertype = computer;
        case 'domain'
            [s,name]=unix('dnsdomainname');
            if uint8(name(end))==10
                name = name(1:end-1);
            end
            machine.domainname = name;
    end
end