function best_subset = ml_read_sasstepdiscout( sasfilename)
% function best_subset = ml_read_sasstepdiscout( sasfilename)
% Parse the output from SAS STEPDISC
% return a list of "best" features

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

% First extract the piece of text containing the table with the
% actual results
NEWLINE = 10;
PAGEBREAK = 12;
f = fopen( sasfilename);
all_text = char(fread( f, [1 inf], 'char'));
summary_pos = findstr( all_text, 'Stepwise Selection Summary');
if length(summary_pos) > 0
    summary = all_text( summary_pos(1):end);
    %clear all_text;
    data = [];
    lastpage = length(summary_pos);
    for page_no = 1:lastpage % For each Summary page do
	last_column_label_pos = findstr( summary, 'ASCC');
	summary = summary( last_column_label_pos(1):end);
	newline_pos = findstr( summary, NEWLINE);
	data_start = newline_pos(2) + 1; % Data starts two lines after
	% the column headings
	if( page_no == lastpage)
	    page = summary( data_start:end);
	    page_end = length(summary);
	else
	    %page_end = findstr( summary, PAGEBREAK);
	    next_page_start = findstr(summary,'The SAS System');
	    current_page_newlines = find( newline_pos < next_page_start(1));
	    page_end = newline_pos( current_page_newlines(end));
	    page = summary( data_start:page_end(1)-1);
	end
	data = [data page];
	summary = summary( page_end(1):end);
    end

    % Then read that table and figure out the best features
    % Write the data to a temporary file
    tempdir = '/tmp';
    [f,filename] = ml_fopentemp( tempdir);
    fullpath = [tempdir '/' filename];
    fwrite(f,data,'char');
    fclose(f);
    % Read the data into a matrix
    [ Step, NumberIn, FeatName, PartialRSquare, FValue, ProbabilityF, ...
      WilksLambda, ProbabilityLambda, AvgSqCanCorr, ProbabilityASCC] ...
	= textread( fullpath, '%d %d %s %f %f %s %f %s %f %s');
    unix(['rm ' fullpath]);
    % Figure out which ones were Entered, and which ones Removed
    RemovedEntered = [1; diff(NumberIn)];  % 1 = Entered, -1 = Removed
    % Turn Feature names into Feature numbers
    for i = 1:length(FeatName)
	FeatNumber(1,i) = str2num(FeatName{i}(4:end))-1;
    end
    % Select the "Good" features as being the ones having Less than 
    % 0.0001 probability of getting an F value of FValue or larger
    % if the null hypthesis (that the classes are the same) was true.
    % I.e. we just go down the list to the place where the first one
    % with ProbabilityF greater than <.0001 is.
    good_feature_idx = [];
    for i=1:length(ProbabilityF)
	PrF = ProbabilityF{i};
	switch( PrF)
	 case '<.0001'
	    good_feature_idx = [good_feature_idx i];
	 otherwise
	end
    end
    best_subset = FeatNumber( good_feature_idx);
    % Finally, sort these "best" features by their F-value or Wilks' lambda
    fval=FValue( good_feature_idx);
    [SortedF,idx]=sort(fval);
    best_subset=best_subset(flipdim(idx,1));
else 
    disp('No features selected by SDA');
    best_subset=[];
end
