%~ Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli  

%~ This file is part of ONSAS.

%~ ONSAS is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation, either version 3 of the License, or
%~ (at your option) any later version.

%~ ONSAS is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%~ GNU General Public License for more details.

%~ You should have received a copy of the GNU General Public License
%~ along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.



% ---------------------------------------------------
% Function for computation of (approximate) deformed configuration of element 
% Inputs:
%  - coordsElem: column vector with the coordinates of the nodes of the element at the reference configuration
%  - dispsElem: colum vector with the displacements of the nodes of the element (rotations are measured with respect to the initial configuration)
%  - elemType: type of element 
%  - locglomat: matrix for transformation from local to global systems
% the element is assumed to be straigt in the reference configuration
%
% Outputs:
%  - vector with xs, ys and zs coordinates for plot.
% ---------------------------------------------------

function [xsdef, ysdef, zsdef, titax, titay, titaz, conecElem ] = outputFrameElementPlotLin( coordsElem, dispsElem, elemType, locglomat )

% reads input
xsref = coordsElem( [ 1   7  ] ) ;
ysref = coordsElem( [ 1+2 7+2] ) ;
zsref = coordsElem( [ 1+4 7+4] ) ;

ndofpnode =  6 ;
conecElem = [] ;

% converts global coord displacements to local
if elemType == 1

  xsdef = xsref + dispsElem( [ 1   1+ndofpnode   ] ) ;
  ysdef = ysref + dispsElem( [ 1+2 1+ndofpnode+2 ] ) ;
  zsdef = zsref + dispsElem( [ 1+4 1+ndofpnode+4 ] ) ;
  titax = [] ;
  conecElem = [1 2] ;

elseif elemType == 2

  nPlotPoints    = 30 ; 

  R = RotationMatrix ( ndofpnode, locglomat ) ;

  localUelem = R' * dispsElem ;

  l      = sqrt( sum( (xsref(2)-xsref(1))^2 + (ysref(2)-ysref(1))^2  + (zsref(2)-zsref(1))^2 ) ) ;
  
  xsglo  = linspace( xsref(1) , xsref(2), nPlotPoints ) ;
  ysglo  = linspace( ysref(1) , ysref(2), nPlotPoints ) ;
  zsglo  = linspace( zsref(1) , zsref(2), nPlotPoints ) ;
  
  xsloc  = linspace( 0 , l, nPlotPoints )' ;
  
  ux     = zeros( size(xsloc) ) ; 
  uy     = zeros( size(xsloc) ) ; 
  uz     = zeros( size(xsloc) ) ;
  titax	 = zeros( size(xsloc) ) ;
  titay	 = zeros( size(xsloc) ) ;
  titaz	 = zeros( size(xsloc) ) ;
  
  LocAxialdofs  = [ 1  1+ndofpnode                        ] ;
  LocTorsidofs	=	[	2	 2+ndofpnode												]	;
  LocBendXYdofs = [ 3  6           3+ndofpnode 6+ndofpnode] ;
  LocBendXZdofs = [ 5  4           5+ndofpnode 4+ndofpnode] ;
  
  for i=1:nPlotPoints,
  
    % bending
    N = bendingInterFuns( xsloc(i), l, 0 );
    uy(i)  = N * localUelem( LocBendXYdofs ) ;

    Rb = eye(4);
    Rb(2,2) = -1;
    Rb(4,4) = -1;

    uz(i)  = N * Rb' * localUelem( LocBendXZdofs ) ;
    
    % stretching
    Nlin1 = (l - xsloc(i))/l ;
    Nlin2 = xsloc(i)/l       ;
    ux(i) = Nlin1 * localUelem(LocAxialdofs(1)) + Nlin2 * localUelem(LocAxialdofs(2)) ;
		
		% torsion
		titax(i) = Nlin1 * localUelem(LocTorsidofs(1)) + Nlin2 * localUelem(LocTorsidofs(2)) ;
		
    % rots
    Nrots = bendingInterFuns( xsloc(i), l, 1 ) ;
    titaz(i) = Nrots * localUelem( LocBendXYdofs ) ; 
    titay(i) = Nrots * Rb' * localUelem( LocBendXZdofs ) ; 
    
		if i<nPlotPoints
			conecElem = [conecElem ; (i-1)+1 (i-1)+2] ;
		end
  end
  
  XsLoc = [ xsloc+ux uy uz];
  
  XsGlo = locglomat * XsLoc' ;
  
  xsdef = ( xsref(1) + XsGlo(1,:) )' ;
  ysdef = ( ysref(1) + XsGlo(2,:) )' ;
  zsdef = ( zsref(1) + XsGlo(3,:) )' ;

end



