% Copyright (C) 2021, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
%   Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini, Sebastian Toro
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

function [ iniNodes, midNodes, endNodes, secc ] = crossSectionVtkSolidConnec( elemCrossSecParams )

  if elemCrossSecParams(1) == 1  || elemCrossSecParams(1) == 2 % general or rectangular section
    if elemCrossSecParams(1) == 1
      % equivalent square section using A = wy * wz
      auxh = sqrt( elemCrossSecParams(2) ) ;   auxb = auxh ;
      secc = [ 12 auxb auxh ] ;
    else
      secc = [ 12 elemCrossSecParams(2) elemCrossSecParams(3) ] ;
    end

		iniNodes = [ 1 2 3 4  ] ;
		midNodes = [          ] ;
    endNodes = [ 5 6 7 8  ] ;

  elseif elemCrossSecParams(1) == 3 % circular section

    secc = [ 25 elemCrossSecParams(2) ] ;
		iniNodes = [ 1 2 3 4 9 10 11 12  ] ;
		midNodes = [ 17 18 19 20         ] ;
    endNodes = [ 5 6 7 8 13 14 15 16 ] ;
  end
