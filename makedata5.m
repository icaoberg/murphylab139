% This script is for cell level classification with purality voting scheme
% The whole computation will takes for about 40 days if only running them on one CPU.   
% Running on clusters is strongly recommended.

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

seednum = 25;
% Part 5-1, corresponding script: "sh makedata5_createnfoldscriptFlag3.sh"
% Usage: Create index file: it takes about one hour for each seed, which this should be done first.
for seedidx = 1:seednum  
   createnfoldscriptR(3,20,seedidx,1,140);
end


% Part 5-2, corresponding script: "sh makedata5_classifyscriptFlag3R.sh"
% Usage: Do classification on 20 class: it takes about 7*6 hours per seed, which can run after index files have been created.
for seedidx = 1:seednum
    classifyscriptR(3,20,seedidx,1,140,1,1,11,5,0,1);   
end


% Part 5-3, corresponding script: "sh makedata5_combinefoldcriptFlag3.sh"
% Usage: Combine results from each fold per seed: it does not take long, which can run after classification results are available.
for seedidx = 1:seednum
    combinescript(3,20,seedidx,1,140,1,1,11,5,0,1);
end


% Part 5-4, corresponding script: "sh makedata5_combineseedcriptFlag3.sh"
% Usage: Combine results from each seed: it takes about 10 minutes, which can run after each classification fold is combined.
combineseedscript(3,20,seednum,1,140,1,1,11,5,0,1);


% Part 5-5, corresponding script: "sh makedata5_combinefoldcriptFlag41.sh"
% Usage: Do classification on ambiguous class: it takes about two hour per seed, which can run after classifiers are trained.
for seedidx = 1:seednum
    classifyotherscript(3,20,seedidx,1,140,1,1,11,5,0,1,41);
end


% Part 5-6, corresponding script: "sh makedata5_combinefoldcriptFlag42.sh"
% Usage: Do classification on punctate_composite class: it takes about one hour per seed, which can run after classifiers are trained.
for seedidx = 1:seednum
    classifyotherscript(3,20,seedidx,1,140,1,1,11,5,0,1,42);
end


% Part 5-7, corresponding script: "sh makedata5_combineseedcriptFlag41.sh"
% Usage: Combine results from each fold and each seed: it takes about 10 minutes, which can run after each classification fold is combined.
combineseedscript(41,20,25,1,140,1,1,11,5,0,1);

  
% Part 5-8, corresponding script: "sh makedata5_combineseedcriptFlag42.sh"
% Usage: Combine results from each fold and each seed: it takes about 10 minutes, which can run after each classification fold is combined.
combineseedscript(42,20,25,1,140,1,1,11,5,0,1);
