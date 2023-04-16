% SPDX-License-Identifier: MIT
% Copyright (c) 2006 Ross K. Snider.  All rights reserved.
%--------------------------------------------------------------------------
% Description:  isoDalton_element_sym2num.m 
%               returns atomic numbers of the element symbols
%--------------------------------------------------------------------------
% Input:  String of elements, example: element_symbols = 'H He Li Be B C' 
%         and array names = isoDalton_element_symbols_read();
%--------------------------------------------------------------------------
% Output:  	returns atomic numbers of the symbols
%
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
function atomic_numbers = isoDalton_element_sym2num(element_symbols,names)

symbol_index = 1;
[s, element_symbols] = strtok(element_symbols);
symbols{symbol_index} = s;
symbol_index = symbol_index + 1;
while length(symbols) > 0
    [s, element_symbols] = strtok(element_symbols);
    if length(s) > 0
        symbols{symbol_index} = s;
        symbol_index = symbol_index + 1;
    else
        break;
    end
end
symbol_count = symbol_index - 1;
name_count = 112;

atomic_numbers = zeros(1,symbol_count);
for i=1:symbol_count
    found_flag = 0;
    for j=1:name_count
        if strcmp(lower(symbols{i}),lower(names{j}.symbol)) == 1
            atomic_numbers(i) = j;
            found_flag = 1;
            break;
        end
    end
    if found_flag == 0
        error([symbols{i} ' is not an element'])
    end
end    
   
