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

function [ Conecvtk, Nodesvtk, vtkDispMat ] = vtkHexa( Nodes, Conec, sectPar , UG, locglos, dispMat ) ;
  
  by = sectPar(2);    bz = sectPar(3);
  Nodesvtk = [] ; Conecvtk = [] ; vtkDispMat = [] ;
  
  nelems = size(Conec,1);   nnodes = size(Nodes,1) ;
  
  NodesDef  = Nodes + [ UG(1:6:end) UG(3:6:end) UG(5:6:end) ] ;
  rotsMat   =         [ UG(2:6:end) UG(4:6:end) UG(6:6:end) ] ;
  
  for i=1:nelems
		nodeselem = Conec(i,1:2)' ;

		[~, locglos] = beamParameters( NodesDef(nodeselem,:) ) ;
		
    ex = locglos(:,1)' ;
    ey = locglos(:,2)' ;
    ez = locglos(:,3)' ;

    nodes = Conec( i,1:2) ;

    matsec = [ -ey*by*.5-ez*bz*.5 ; ...
               +ey*by*.5-ez*bz*.5 ; ...
               +ey*by*.5+ez*bz*.5 ; ...
               -ey*by*.5+ez*bz*.5 ] ; 
    
    vecrot = rotsMat( nodes(1), :) ;
    matsecR = ( expon( vecrot) * matsec' )' ;

    candsini = NodesDef( nodes(1), : ) + matsecR ;
    

    vecrot = rotsMat( nodes(2), :) ;
    matsecR = ( expon( vecrot) * matsec' )' ;

    candsfin = NodesDef( nodes(2), : ) + matsecR ;
          
    Nodesvtk = [ Nodesvtk ; candsini ; candsfin ] ; 
    Conecvtk = [ Conecvtk ; 5  (i-1)*8+(1:8) ] ;

    vtkDispMat = [ vtkDispMat ; ones(4,1)*dispMat(nodes(1),:) ] ;
    vtkDispMat = [ vtkDispMat ; ones(4,1)*dispMat(nodes(2),:) ] ;
      
  end


end

