% Copyright (C) 2019, Jorge M. Perez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquin Viera, Mauricio Vanzulli  
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

% function that creates the hexahedron element to visualize results with vtk file

function [ Nodesvtk, Conecvtk ] = vtkBeam2SolidConverter( nodesCoords, sectPar, Ue, Rr ) ;

typeSolid = sectPar(1) ;

vecrotNode1 = Ue( 2:2:6  ) ;
vecrotNode2 = Ue( 8:2:12 ) ;

if typeSolid == 12
  by = sectPar(2) ;  bz = sectPar(3) ;
  
  Nodesvtk = zeros( 8,3 ) ;  Conecvtk = zeros( 8,1 ) ;

  [~, locglos] = beamParameters( nodesCoords ) ;
  
  ex = locglos(:,1)' ;
  ey = locglos(:,2)' ;
  ez = locglos(:,3)' ;

  % matrix with coords of four vertices of cross section to be plotted
  matsec = [ -ey*by*.5-ez*bz*.5 ; ...
             +ey*by*.5-ez*bz*.5 ; ...
             +ey*by*.5+ez*bz*.5 ; ...
             -ey*by*.5+ez*bz*.5 ] ;
  
  % rotated section
  matsecR = ( Rr * locglos' * expon( vecrotNode1 ) * matsec' )' ;
  
  candsini = (nodesCoords( 1, : )+Ue(1:2:5)') + matsecR      ;
  
  matsecR = ( Rr * locglos' * expon( vecrotNode2) * matsec' )'  ;

  candsfin = (nodesCoords( 2, : )+Ue(7:2:11)') + matsecR      ;
        
  Nodesvtk = [ candsini ; candsfin ] ; 
  Conecvtk = [ 12 1:8 ] ;    
end

