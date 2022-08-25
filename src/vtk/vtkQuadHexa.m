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
 
% function that creates the quadratic hexahedron element to visualize results with vtk file

function [ Conecvtk, Nodesvtk, vtkDispMat ] = vtkQuadHexa( Nodes, Conec, sectPar , UG, locglos, dispMat ) ;

  R = sectPar(2) ;
  Nodesvtk = [] ; Conecvtk = [] ; vtkDispMat = [] ;

  nelems = size(Conec,1);   nnodes = size(Nodes,1) ;

  NodesDef  = Nodes + [ UG(1:6:end) UG(3:6:end) UG(5:6:end) ] ;
  rotsMat   =         [ UG(2:6:end) UG(4:6:end) UG(6:6:end) ] ;

  for i = 1:nelems
    nodeselem = Conec(i,1:2)' ;
		[~, locglos] = beamParameters( NodesDef(nodeselem,:) ) ;

    ex = locglos(:,1)' ;
    ey = locglos(:,2)' ;
    ez = locglos(:,3)' ;

    nodeselem = Conec(i,1:2) ;
    %
    matSecNodesExtVertices   = [ -ey*R/sqrt(2)-ez*R/sqrt(2) ; ...
                                 +ey*R/sqrt(2)-ez*R/sqrt(2) ; ...
                                 +ey*R/sqrt(2)+ez*R/sqrt(2) ; ...
                                 -ey*R/sqrt(2)+ez*R/sqrt(2) ] ;

    matSecNodesExtInter      = [            0-ez*R       ; ...
                                       +ey*R-0           ; ...
                                            0+ez*R       ; ...
                                       -ey*R+0           ] ;
    % Definida solo por claridad ...
    matSecNodesInter         = [ -ey*R/sqrt(2)-ez*R/sqrt(2) ; ...
                                 +ey*R/sqrt(2)-ez*R/sqrt(2) ; ...
                                 +ey*R/sqrt(2)+ez*R/sqrt(2) ; ...
                                 -ey*R/sqrt(2)+ez*R/sqrt(2) ] ;

    % Rot section 1
    vecrot = rotsMat( nodeselem(1), :) ;
    matSecNodesExtVerticesR  = ( expon( vecrot) * matSecNodesExtVertices' )'  ;
    matSecNodesExtInterR     = ( expon( vecrot) * matSecNodesExtInter' )'     ;

    nodeExtVertices1  = NodesDef(nodeselem(1),:) + matSecNodesExtVerticesR    ;
    nodeExtInter1     = NodesDef(nodeselem(1),:) + matSecNodesExtInterR       ;

    % Rot section 2
    vecrot = rotsMat( nodeselem(2), :) ;
    matSecNodesExtVerticesR  = ( expon( vecrot) * matSecNodesExtVertices' )'  ;
    matSecNodesExtInterR     = ( expon( vecrot) * matSecNodesExtInter' )'     ;

    nodeExtVertices2  = NodesDef(nodeselem(2),:) + matSecNodesExtVerticesR     ;
    nodeExtInter2     = NodesDef(nodeselem(2),:) + matSecNodesExtInterR        ;

    % Rot intermediate section
    vecrot = ( rotsMat( nodeselem(1), :) + rotsMat( nodeselem(2), :) ) / 2            ;
    matSecNodesInterR  = ( expon( vecrot) * matSecNodesInter' )'               ;
    nodeInt           = ( NodesDef(nodeselem(1),:) + NodesDef(nodeselem(2),:) ) / 2 + matSecNodesInterR ; % Interpolated linearly ...

    % Nodes
    Nodesvtk = [ Nodesvtk ; nodeExtVertices1 ; nodeExtVertices2 ; nodeExtInter1 ; nodeExtInter2 ; nodeInt ] ;

    % Connectivity
    Conecvtk = [ Conecvtk ; 25 (i-1)*20+(1:20) ] ;

    % Disps mat
    vtkDispMat = [ vtkDispMat ; ones(4,1)*dispMat(nodeselem(1),:) ] ;
    vtkDispMat = [ vtkDispMat ; ones(4,1)*dispMat(nodeselem(2),:) ] ;
    vtkDispMat = [ vtkDispMat ; ones(4,1)*dispMat(nodeselem(1),:) ] ;
    vtkDispMat = [ vtkDispMat ; ones(4,1)*dispMat(nodeselem(2),:) ] ;
    vtkDispMat = [ vtkDispMat ; ones(4,1)*(dispMat(nodeselem(1),:)+(dispMat(nodeselem(1),:)))/2 ] ;

  end

end
