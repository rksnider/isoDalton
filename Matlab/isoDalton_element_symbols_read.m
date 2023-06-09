% SPDX-License-Identifier: MIT
% Copyright (c) 2006 Ross K. Snider.  All rights reserved.
%--------------------------------------------------------------------------
% Description:  isoDalton_element_symbols_read.m 
%               This file parses the atomic name and symbol file 
%               isoDalton_element_symbols.txt  The data file was created 
%               from http://www.chem.qmul.ac.uk/iupac/AtWt/
%--------------------------------------------------------------------------
% Input:  None (reads a file that should be present)
%--------------------------------------------------------------------------
% Output:  	The cell array element with the following fields:
%               element{atomic_number}.symbol    % atomic_symbol;
%               element{atomic_number}.name      % element name
%--------------------------------------------------------------------------
% This software is associated with the following paper:
% Snider, R.K. Efficient Calculation of Exact Mass Isotopic Distributions
% J Am Soc Mass Spectrom 2007, Vol 18/8 pp. 1511-1515.
% The digital object identifier (DOI) link to paper:  
% http://dx.doi.org/10.1016/j.jasms.2007.05.016
%--------------------------------------------------------------------------
% Author:       Ross K. Snider
% Company:      Montana State University
% Create Date:  April 27, 2006
% Revision:     1.0
% License: MIT  (opensource.org/licenses/MIT)
%--------------------------------------------------------------------------
function names = isoDalton_element_symbols_read()

fid=fopen('isoDalton_element_symbols.txt');
if fid < 0
    error('Cannot open isoDalton_element_symbols.txt');
end

line = fgetl(fid);

skip_header = 1;    % skip header  (lines starting with %)
while skip_header == 1 & length(line) > 0
    if line(1) == '%';
        line = fgetl(fid);
    else
        skip_header = 0;
    end
end

while length(line) == 0   % skip blank lines
    line = fgetl(fid);
end


skip_flag = 1;
while 1
    if skip_flag == 0
        line = fgetl(fid);
        if ~ischar(line), break, end
    end

    [w, line] = strtok(line);
    [s, line] = strtok(line);
    [n, line] = strtok(line);
   
   i = str2num(w);
   names{i}.symbol = s;
   names{i}.name = n;

   skip_flag = 0;    
end

fclose(fid);



