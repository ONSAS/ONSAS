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

% this function creates the conectivities and coordinates data structures for producing
% vtk files of the solids or structures.

function [ Nodesvtk, vtkDispMat, vtkNormalForces, Conecvtk, elem2VTKCellMap ] = ...
  vtkGeometry ( coordsElemsMat, Conec, secc, U, normalForces, Nodes, elementsParamsMat )
  %~ vtkGeometry ( Nodes, Conec, secc, U )
% ------------------------------------------------------------------------------

% filter non-material elements
inds = find( Conec( :, 4+1 )== 0) ;

Conec( inds,:) = [] ;

types = elementsParamsMat( Conec(:, 4+2), 1) ;

nElemsTrussOrFrame = sum( ( types == 1 ) + ( types == 2 ) ) ;
nElemsTetraOrPlate = sum( ( types == 3 ) + ( types == 4 ) ) ;

nelems = size( Conec, 1 ) ;

vtkDispMat = [] ;

if nElemsTrussOrFrame == nelems, % all are truss or beams

  Nodesvtk = [] ;
  vtkNormalForces = [] ;
  
  if nargout>1

    counterNodesVtk = 0 ;
    counterCellsVtk = 0 ;
    Conecvtk        = [] ;

    % mapping from the original connectivity to the vtk cells
    elem2VTKCellMap = (1:nelems)' ; % default: i-to-i . Colums should be changed.
  end
  
  % loop in elements
  for i = 1:nelems
    elemMat  = Conec(i, 4+1 ) ;
    elemType = elementsParamsMat( Conec(i, 4+2 ), 1 ) ;
    elemSec  = Conec(i, 4+5 ) ;
        
    nodeselem  = Conec(i,1:2)' ;
    dofselem   = nodes2dofs( nodeselem , 6 ) ;
    dispsElem  = U( dofselem ) ;
      
    coordsElem         = coordsElemsMat( i, 1 : 12 )'  ;
    %~ coordsElem(1:2:6 ) = Nodes( nodeselem(1), : )' ;      
    %~ coordsElem(7:2:12) = Nodes( nodeselem(2), : )' ;
    
    %~ dispsElem     = zeros(12,1) ;
    %~ dispsElem( 6) =  .5 ;
    %~ dispsElem(12) = -.5 ;
    
    [ xdef, ydef, zdef, conecElem, titax, titay, titaz, Rr ] = outputFrameElementPlot( coordsElem, dispsElem, elemType ) ;
 
    % nodes sub-division 
    ndivNodes = conecElem( end, 2 ) ;
      
    % Local x coords in the reference config
    xloc = linspace( coordsElem(1), coordsElem( 7), ndivNodes )' ;  
    yloc = linspace( coordsElem(3), coordsElem( 9), ndivNodes )' ;  
    zloc = linspace( coordsElem(5), coordsElem(11), ndivNodes )' ;

    if length( secc ) > 0
        

      for j = 1 : (ndivNodes-1)

        dispsSubElem       = zeros(12,1) ;
        dispsSubElem(1:6 ) = [ xdef(j  )-xloc(j  ) titax(j  ) ydef(j  )-yloc(j  ) titay(j  ) zdef(j  )-zloc(j  ) titaz(j  ) ]' ;
        dispsSubElem(7:12) = [ xdef(j+1)-xloc(j+1) titax(j+1) ydef(j+1)-yloc(j+1) titay(j+1) zdef(j+1)-zloc(j+1) titaz(j+1) ]' ;

        coordSubElem             = [ xloc(j:(j+1))  yloc(j:(j+1))  zloc(j:(j+1)) ] ;
        
        [ NodesCell, ConecCell ] = vtkBeam2SolidConverter ( coordSubElem, secc, dispsSubElem, Rr ) ;
        
        if j==1
          Nodesvtk = [ Nodesvtk ; NodesCell( 1:4, :) ] ;
          vtkDispMat = [ vtkDispMat ; ones(4,1)*(dispsSubElem(1:2:5)') ] ;
        end
        
        Nodesvtk = [ Nodesvtk ; NodesCell( 5:end, : ) ] ; 
        vtkDispMat = [ vtkDispMat ; ones(4,1)*(dispsSubElem(7:2:11)') ] ;

        if nargout > 3
          ConecCell(:,2:end) = ConecCell(:,2:end) + (j-1) * 4 + counterNodesVtk ;
          Conecvtk = [ Conecvtk ; ConecCell ] ;
        end

      end

      vtkNormalForces = [ vtkNormalForces ; ones( (ndivNodes-1),1)*normalForces(i) ] ;


      if nargout > 3
        elem2VTKCellMap( i, 1:9 ) = (1:(ndivNodes-1)) + counterCellsVtk ;
      
        counterNodesVtk = counterNodesVtk +   ndivNodes * 4   ;
        counterCellsVtk = counterCellsVtk + ( ndivNodes - 1 ) ;
      end

    else
      %~ polilyine
      
      %~ % Undeformed nodes matrix
      %~ NodesAux 		= [ NodesAux			; xloc yloc zloc ] ;   
  
      %~ % Deformed nodes matrix
      %~ NodesDefAux = [ NodesDefAux		; xdef ydef zdef ] ;
  
      %~ % Rots aux
      %~ RotsAux     = [ RotsAux       ; titax titay titaz ] ;   

      %~ % New connectivity
      %~ if elemType == 2
        %~ ConecAux = [ ConecAux ; ((i-1)*ndivselem+conecElem) zeros(size(conecElem,1),2) ...
                                %~ ones(size(conecElem,1),1)*elemMat ones(size(conecElem,1),1)*elemSec ones(size(conecElem,1),1)*elemType ] ; 
  
      %~ else 
        %~ ConecAux = [ ConecAux ; 2*i-1 2*i zeros(size(conecElem,1),2) ...
                                %~ ones(size(conecElem,1),1)*elemMat ones(size(conecElem,1),1)*elemSec ones(size(conecElem,1),1)*elemType ] ;  
      %~ end
        
    end
  end
  
elseif nElemsTetraOrPlate == nelems  % all are tetraedra or plates

  vtkNormalForces = [] ;
  vtkDispMat = reshape( U(1:2:end)', [3, size(Nodes,1) ])' ;

  Nodesvtk = Nodes + vtkDispMat ;
    
  if nargout > 3
    Conecvtk = [ types Conec(:,1:4) ] ;
    elem2VTKCellMap = (1:nelems)' ; % default: i-to-i . Colums should be changed.
  end

else
  error(' mixed elements plot not implemented yet. Please create an issue.')

end
