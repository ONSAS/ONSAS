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

checkFields(materials, expectedFields)

end
function checkFields(mystruct, expectedFields)

    myfields = fieldnames(mystruct)
    len = length(chars) ; % as cell
    for i = 1:len
        chars{i}
        fields{i}
        if strcmp(chars{i},fields{i}) ~= 1
            error("The struct %s does not contain the struct %s, to see the fields please check: https://onsas.github.io/ONSAS.m/dev/howtouse/creatingModels/", structName, fields{i}) 
        end     
    end

end