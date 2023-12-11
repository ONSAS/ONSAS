% Copyright 2023, ONSAS Authors (see documentation)
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


function elements = setDefauAeroCrossSecProps( elements )

  for i=1:length(elements)
    if iscell(elements(i).elemCrossSecParams)
      cross_sec_name = elements(i).elemCrossSecParams{1} ;
      
      if strcmp( cross_sec_name, 'circle' ) || strcmp( cross_sec_name, 'pipe' )
        if isempty( elements(i).chordVector)
          elements(i).chordVector = [ 0 0 elements(i).elemCrossSecParams{2} ] ;
        end
      
        if iscell( elements(i).aeroCoefFunctions)
          if isempty( elements(i).aeroCoefFunctions{1}) && isempty( elements(i).aeroCoefFunctions{2}) && isempty( elements(i).aeroCoefFunctions{3})
            anonymus_null = @(beta,Re) 0 ;
            elements(i).aeroCoefFunctions = { 'innerDragCoefCircular'; anonymus_null; anonymus_null } ; 
          end
        end

      end
    end
  end
