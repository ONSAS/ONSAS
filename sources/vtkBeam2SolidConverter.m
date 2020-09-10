% Copyright (C) 2020, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
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

% function that creates the hexahedron element to visualize results with vtk file

function [ Nodesvtk, Conecvtk ] = vtkBeam2SolidConverter( nodesCoords, sectPar, Ue, Rr ) ;

typeSolid = sectPar(1) ;

vecrotNode1 = Ue( 2:2:6  ) ;
vecrotNode2 = Ue( 8:2:12 ) ;

if typeSolid == 12 % vtkHexa

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

  candsini = ones(4,1) * (nodesCoords( 1, : )+Ue(1:2:5)') + matsecR      ;
  
  matsecR = ( Rr * locglos' * expon( vecrotNode2) * matsec' )'  ;

  candsfin = ones(4,1) * (nodesCoords( 2, : )+Ue(7:2:11)') + matsecR      ;
        
  Nodesvtk = [ candsini ; candsfin ] ; 
  Conecvtk = [ 12 1:8 ] ; 
  
elseif typeSolid == 25 % vtkQuadHexa

	R = sectPar(2) / 2 ;
  Nodesvtk = [] ; Conecvtk = [] ;
  
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
	
	nodeExtVertices1  = NodesDef(1,:) + matSecNodesExtVerticesR    ;
	nodeExtInter1     = NodesDef(1,:) + matSecNodesExtInterR       ;   
    
	% Rot section 2
	matSecNodesExtVerticesR  = ( expon( vecrotNode2 ) * matSecNodesExtVertices' )'  ;
	matSecNodesExtInterR     = ( expon( vecrotNode2 ) * matSecNodesExtInter' )'     ;
	
	nodeExtVertices2  = NodesDef(2,:) + matSecNodesExtVerticesR     ;
	nodeExtInter2     = NodesDef(2,:) + matSecNodesExtInterR        ;  
	
  % Rot intermediate section
	vecrot = ( vecrotNode1 + vecrotNode2 ) / 2            ;
	matSecNodesInterR = ( expon( vecrot) * matSecNodesInter' )'               ;
	nodeInt           = ( NodesDef(1,:) + NodesDef(2,:) ) / 2 + matSecNodesInterR ; % Interpolated linearly ...
	
  % Nodes
	Nodesvtk = [ Nodesvtk         ; ...
               nodeExtVertices1 ; ...
               nodeExtVertices2 ; ...
               nodeExtInter1    ; ...
               nodeExtInter2    ; ...
               nodeInt          ] ;

	% Connectivity
	Conecvtk = [ Conecvtk ; 25 1:20 ] ;  
     
end

