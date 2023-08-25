% Copyright 2023, Joaquin Viera, Jorge M. Perez Zerpa.
%
% This file is part of ONSAS.
%
% ONSAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ONSAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

function checkONSASFields( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams )

    checkFields(materials, {'hyperElasModel', 'hyperElasParams','density'})

function checkFields(mystruct, expectedFields)

    myfields = fieldnames(mystruct)
    len = length(expectedFields) ; % as cell
    
    for i = 1:length(myfields)
        j = 1;
        while (j <= len) && (strcmp(myfields{i},expectedFields{j}) ~= 1)
            j = j+1;
        end
        if j == (len+1)
            error('The field %s is not a correct field name', myfields{i}) 
        end     
    end
