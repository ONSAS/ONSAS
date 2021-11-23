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

  dispIniSection = dispsElem(1:2:5)  ;
  dispEndSection = dispsElem(7:2:11) ;

  if typeSolid == 12 % vtkHexa

     by = sectPar(2) ;  bz = sectPar(3) ;

     Nodesvtk = zeros( 8,3 ) ;  Conecvtk = zeros( 8,1 ) ;

     % coordinates of nodes 1 and 2 in deformed configuration
%5     node1Def = nodesCoords(1:3) + nodalDisp(1:3) ;
%     node2Def = nodesCoords(4:6) + nodalDisp(4:6) ;

%     locglos = beamRefConfRotMat( (node2Def - node1Def ) ) ;

%     ex = locglos(:,1)' ;
%     ey = locglos(:,2)' ;%
%     ez = locglos(:,3)' ;
     ex = [1 0 0] ;
     ey = [0 1 0] ;
     ez = [0 0 1] ;

    % 1- compute the vectors of the section in Rr coords
    % matrix with coords of four vertices of cross section to be plotted
    matrixSectionIni = [ -ey*by*.5-ez*bz*.5 ; ...
                         +ey*by*.5-ez*bz*.5 ; ...
                         +ey*by*.5+ez*bz*.5 ; ...
                         -ey*by*.5+ez*bz*.5 ] ;
    matrixSectionEnd = matrixSectionIni ;
    matrixSectionMed = [] ;
  end

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

  NodesRefvtk = [ matrixRefIni ; matrixRefEnd ] ;

  Dispsvtk = Nodesvtk - NodesRefvtk ;


return
	R = sectPar(2) / 2 ;
  NodesDef  = nodesCoords + [ Ue(1:6:end) Ue(3:6:end) Ue(5:6:end) ] ;
  rotsMat   = 			        [ Ue(2:6:end) Ue(4:6:end) Ue(6:6:end) ] ;

  [~, locglos] = beamParameters( NodesDef ) ;


  ex = locglos(:,1)' ;
	ey = locglos(:,2)' ;
	ez = locglos(:,3)' ;

	%
	matSecNodesExtVertices   = [ -ey*R/sqrt(2)-ez*R/sqrt(2) ; ...
															 +ey*R/sqrt(2)-ez*R/sqrt(2) ; ...
															 +ey*R/sqrt(2)+ez*R/sqrt(2) ; ...
															 -ey*R/sqrt(2)+ez*R/sqrt(2) ] ;

	matSecNodesInter         = matSecNodesExtVertices ;

	matSecNodesExtInter      = [             0-ez*R       ; ...
                                       +ey*R-0           ; ...
																				   0+ez*R       ; ...
                                       -ey*R+0           ] ;


	% Rot section 1
	matSecNodesExtVerticesR  = ( expon( vecrotNode1 ) * matSecNodesExtVertices' )'  ;
	matSecNodesExtInterR     = ( expon( vecrotNode1 ) * matSecNodesExtInter' )'     ;

	nodeExtVertices1  = ones(4,1) * NodesDef(1,:) + matSecNodesExtVerticesR    ;
	nodeExtInter1     = ones(4,1) * NodesDef(1,:) + matSecNodesExtInterR       ;

	% Rot section 2
	matSecNodesExtVerticesR  = ( expon( vecrotNode2 ) * matSecNodesExtVertices' )'  ;
	matSecNodesExtInterR     = ( expon( vecrotNode2 ) * matSecNodesExtInter' )'     ;

	nodeExtVertices2  = ones(4,1) * NodesDef(2,:) + matSecNodesExtVerticesR     ;
	nodeExtInter2     = ones(4,1) * NodesDef(2,:) + matSecNodesExtInterR        ;

  % Rot intermediate section
	vecrot = ( vecrotNode1 + vecrotNode2 ) / 2            ;
	matSecNodesInterR = ( expon( vecrot) * matSecNodesInter' )'               ;
	nodeInt           = ones(4,1) * ( NodesDef(1,:) + NodesDef(2,:) ) / 2 + matSecNodesInterR ; % Interpolated linearly ...

  % Nodes
	Nodesvtk = [ Nodesvtk         ; ...
               nodeExtVertices1 ; ...
               nodeExtVertices2 ; ...
               nodeExtInter1    ; ...
               nodeExtInter2    ; ...
               nodeInt          ] ;

	% Connectivity
	Conecvtk = [ Conecvtk ; 25 1:20 ] ;

%end
