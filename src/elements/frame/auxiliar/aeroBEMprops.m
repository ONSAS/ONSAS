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

function [ chordVector, aoa, cl, cd, cm, fstat, clinv, clfullysep ] = aeroBEMprops( aeroAirfoils, sectionData ) ;

for i = length(aeroAirfoils)
    [aoa(:,i), cl(:,i), cd(:,i), cm(:,i), fstat(:,i), clinv(:,i), clfullysep(:,i)] = readairfoildata( aeroAirfoils(i) );
end

twist       = sectionData(:,2); % Angle in DEGRESS
chord       = sectionData(:,3);

for j = 1:length(twist)
    chordVector(1,j) = [0, 0 , chord(j)] ;
end
end