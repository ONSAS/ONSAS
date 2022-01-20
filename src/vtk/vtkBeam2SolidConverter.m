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

%md function that creates the hexahedron element to visualize results with vtk file
function [ Nodesvtk, Conecvtk, Dispsvtk ] = vtkBeam2SolidConverter( ...
  coordsElemNodes, dispsElem, coordLocSubElem, dispLocIni, dispLocEnd, ...
  locRotIni, locRotEnd, sectPar, Rr, R0 ) ;

typeSolid = sectPar(1) ;

dispIniSection = dispsElem(1:2:5) ;    dispEndSection = dispsElem(7:2:11) ;

dispLocMed = 0.5 * ( dispLocIni + dispLocEnd ) ;

ex = [1 0 0] ;
ey = [0 1 0] ;
ez = [0 0 1] ;

if typeSolid == 12 % vtkHexa

   by = sectPar(2) ;  bz = sectPar(3) ;

   Nodesvtk = zeros( 8,3 ) ;  Conecvtk = zeros( 8,1 ) ;

  % 1- compute the vectors of the section in Rr coords
  % matrix with coords of four vertices of cross section to be plotted
  matrixSectionIni = [ -ey*by*.5-ez*bz*.5 ; ...
                       +ey*by*.5-ez*bz*.5 ; ...
                       +ey*by*.5+ez*bz*.5 ; ...
                       -ey*by*.5+ez*bz*.5 ] ;
  matrixSectionEnd = matrixSectionIni ;

  % 2- apply the local rotation for the initial and final sections
  matrixRotatedSectionIni = ( expon( locRotIni ) * matrixSectionIni' )' ;
  matrixRotatedSectionEnd = ( expon( locRotEnd ) * matrixSectionEnd' )' ;

  % 3- add local displacement and position
  matrixDisplacedSectionIni = matrixRotatedSectionIni + [ coordLocSubElem(1) 0 0] + dispLocIni' ;
  matrixDisplacedSectionEnd = matrixRotatedSectionEnd + [ coordLocSubElem(2) 0 0] + dispLocEnd' ;

  % 4- apply Rr' change basis matrix
  matrixRotatedSectionIni = ( Rr * matrixDisplacedSectionIni' )' ;
  matrixRotatedSectionEnd = ( Rr * matrixDisplacedSectionEnd' )' ;

  defPosIniSec = coordsElemNodes(1:3) + dispIniSection ;

  % 5- add nodal displacements in e1,e2,e3 system
  nodesVtkSectionIni = matrixRotatedSectionIni + defPosIniSec' ;
  nodesVtkSectionEnd = matrixRotatedSectionEnd + defPosIniSec' ;

  Nodesvtk = [ nodesVtkSectionIni ; nodesVtkSectionEnd ] ;
  Conecvtk = [ 12 0:7 ] ; % in vtk indexation (from 0)

  matrixRefIni = ( R0 * (matrixSectionIni + [ coordLocSubElem(1) 0 0] )')' + coordsElemNodes(1:3)' ;
  matrixRefEnd = ( R0 * (matrixSectionEnd + [ coordLocSubElem(2) 0 0] )')' + coordsElemNodes(1:3)' ;

  NodesRefvtk = [ matrixRefIni; matrixRefEnd ] ;

elseif typeSolid == 25 % vtkQuadHexa
  R = sectPar(2) / 2 ;

  Nodesvtk = zeros( 20,3 ) ;  Conecvtk = zeros( 20,1 ) ;

  matrixSectionIni = [ -ey*R/sqrt(2)-ez*R/sqrt(2) ; ...
											 +ey*R/sqrt(2)-ez*R/sqrt(2) ; ...
											 +ey*R/sqrt(2)+ez*R/sqrt(2) ; ...
											 -ey*R/sqrt(2)+ez*R/sqrt(2) ] ;
  matrixSectionEnd = matrixSectionIni ;

  matrixSectionCurvIni = [ 0-ez*R      ; ...
                      +ey*R-0      ; ...
                      0+ez*R       ; ...
                      -ey*R+0      ] ;

  matrixSectionCurvEnd = matrixSectionCurvIni ;

  matrixSectionMed = matrixSectionIni ;

  % 2- apply the local rotation for the initial and final sections
  matrixRotatedSectionIni     = ( expon( locRotIni )                  * matrixSectionIni'     )' ;
  matrixRotatedSectionEnd     = ( expon( locRotEnd )                  * matrixSectionEnd'     )' ;
  matrixRotatedSectionCurvIni = ( expon( locRotIni )                  * matrixSectionCurvIni' )' ;
  matrixRotatedSectionCurvEnd = ( expon( locRotEnd )                  * matrixSectionCurvEnd' )' ;
  matrixRotatedSectionMed     = ( expon( (locRotIni+locRotEnd)*0.5  ) * matrixSectionMed'     )' ;

  % 3- add local displacement and position
  matrixDisplacedSectionIni     = matrixRotatedSectionIni     + [ coordLocSubElem(1) 0 0] + dispLocIni' ;
  matrixDisplacedSectionEnd     = matrixRotatedSectionEnd     + [ coordLocSubElem(2) 0 0] + dispLocEnd' ;
  matrixDisplacedSectionCurvIni = matrixRotatedSectionCurvIni + [ coordLocSubElem(1) 0 0] + dispLocIni' ;
  matrixDisplacedSectionCurvEnd = matrixRotatedSectionCurvEnd + [ coordLocSubElem(2) 0 0] + dispLocEnd' ;
  matrixDisplacedSectionMed     = matrixRotatedSectionMed     + [ (coordLocSubElem(1)+coordLocSubElem(2))*.5 0 0] + dispLocMed' ;

  % 4- apply Rr' change basis matrix
  matrixRotatedSectionIni     = ( Rr * matrixDisplacedSectionIni'     )' ;
  matrixRotatedSectionEnd     = ( Rr * matrixDisplacedSectionEnd'     )' ;
  matrixRotatedSectionCurvIni = ( Rr * matrixDisplacedSectionCurvIni' )' ;
  matrixRotatedSectionCurvEnd = ( Rr * matrixDisplacedSectionCurvEnd' )' ;
  matrixRotatedSectionMed     = ( Rr * matrixDisplacedSectionMed'     )' ;

  defPosIniSec = coordsElemNodes(1:3) + dispIniSection ;

  % 5- add nodal displacements in e1,e2,e3 system
  nodesVtkSectionIni     = matrixRotatedSectionIni     + defPosIniSec' ;
  nodesVtkSectionEnd     = matrixRotatedSectionEnd     + defPosIniSec' ;
  nodesVtkSectionCurvIni = matrixRotatedSectionCurvIni + defPosIniSec' ;
  nodesVtkSectionCurvEnd = matrixRotatedSectionCurvEnd + defPosIniSec' ;
  nodesVtkSectionMed     = matrixRotatedSectionMed     + defPosIniSec' ;

  Nodesvtk = [ nodesVtkSectionIni ;     nodesVtkSectionEnd; ...
               nodesVtkSectionCurvIni;  nodesVtkSectionCurvEnd;  nodesVtkSectionMed ] ;
  Conecvtk = [ 25 0:19 ] ; % in vtk indexation (from 0)

  matrixRefIni = ( R0 * (matrixSectionIni + [ coordLocSubElem(1) 0 0] )')' + coordsElemNodes(1:3)' ;
  matrixRefEnd = ( R0 * (matrixSectionEnd + [ coordLocSubElem(2) 0 0] )')' + coordsElemNodes(1:3)' ;
  matrixRefCurvIni = ( R0 * (matrixSectionCurvIni + [ coordLocSubElem(1) 0 0] )')' + coordsElemNodes(1:3)' ;
  matrixRefCurvEnd = ( R0 * (matrixSectionCurvEnd + [ coordLocSubElem(2) 0 0] )')' + coordsElemNodes(1:3)' ;
  matrixRefMed = ( R0 * (matrixSectionMed + [ (coordLocSubElem(1)+coordLocSubElem(2))*.5 0 0] )')' + coordsElemNodes(1:3)' ;

  NodesRefvtk = [ matrixRefIni ;  matrixRefEnd; matrixRefCurvIni; matrixRefCurvEnd; matrixRefMed ] ;

end

Dispsvtk = Nodesvtk - NodesRefvtk ;
