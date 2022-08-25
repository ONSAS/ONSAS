% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% Joaquin Viera, J. Bruno Bazzano, Marcelo Forets, Jean-Marc Battini.
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

  crossSecShape  = elemCrossSecParams{1} ;
  crossSecParams = elemCrossSecParams{2} ;
  if size( elemCrossSecParams, 1 ) > 2
    crossSecBoundary = elemCrossSecParams{3} ;
    if ~strcmp( crossSecShape, 'generic' )
      error('check element sections')
    end
  end

  if strcmp( crossSecShape, 'circle' )
    diameter = crossSecParams(1) ;
    secc = [ 25 diameter ] ;
    iniNodes = [ 1 2 3 4 9 10 11 12  ] ;
    midNodes = [ 17 18 19 20         ] ;
    endNodes = [ 5 6 7 8 13 14 15 16 ] ;

  elseif strcmp( crossSecShape, 'rectangle' ) || ( strcmp( crossSecShape, 'generic' ) && ~exist('crossSecBoundary') )

    if strcmp( crossSecShape, 'generic' )
      area = crossSecParams(1) ;
      % equivalent square section using A = wy * wz
      auxh = sqrt( area ) ;   auxb = auxh ;
      secc = [ 12 auxb auxh ] ;
    else
      by = crossSecParams(1);  bz = crossSecParams(2);
      secc = [ 12 by bz ] ;
    end

		iniNodes = [ 1 2 3 4  ] ;
		midNodes = [          ] ;
    endNodes = [ 5 6 7 8  ] ;

  else
    iniNodes = [ ] ;
    midNodes = [          ] ;
    endNodes = [   ] ;

  end
