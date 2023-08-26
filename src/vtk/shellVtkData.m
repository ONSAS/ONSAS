% Copyright 2022, Jorge M. Perez Zerpa.
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
 
function [ vtkNodes, vtkConec, vtkNodalDisps ] ...
   = shellVtkData( Nodes, Conec, elemCrossSecParams, U )

  vtkNodes        = [] ;
  vtkConec        = [] ;
  vtkNodalDisps   = [] ;

  nPlotSubElements = 10 ; % number of plot subsegments
  counterNodes     = 0 ;
  nelem            = size(Conec,1) ;

  % thickness
  tz = elemCrossSecParams{2};


  for i=1:nelem

    % nodes and degrees of freedom of current element
    nodesElem  = Conec(i,1:3)               ;
    dofsElem   = nodes2dofs( nodesElem, 6 ) ;

    crossVector = cross( ...
        Nodes( nodesElem(2),:) - Nodes( nodesElem(1),:) , ...
        Nodes( nodesElem(3),:) - Nodes( nodesElem(1),:) ...
        ) ;
    normalVector = crossVector / norm(crossVector) ;

    Nodesvtk =   [ Nodes( nodesElem,:) ; ...
                   Nodes( nodesElem,:) ] ... 
               + tz*.5*[ ones(3,1);-ones(3,1)] * normalVector ;

    % column vector with displacements of the dofs of the current element
    dispsElem  = U( dofsElem ) ;

    aux = reshape( dispsElem(1:2:end)', [3, 3 ])' ;

    Dispsvtk = [aux; aux];

    Conecvtk = [ 13  (nodes2dofs(i,6)-1)' ] ;

    vtkNodes             = [ vtkNodes ;     Nodesvtk ] ;
    vtkConec             = [ vtkConec ;     Conecvtk ] ;
    vtkNodalDisps        = [ vtkNodalDisps; Dispsvtk ] ;

  end % for elements

