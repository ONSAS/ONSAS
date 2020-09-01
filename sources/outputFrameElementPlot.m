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


% ---------------------------------------------------
% Function for computation of (approximate) deformed configuration of one 
% beam or truss element using interpolations.
% 
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

function [ xsdef, ysdef, zsdef, conecElem, titax, titay, titaz, Rr ] = ...
  outputFrameElementPlot( coordsElem, dispsElem, elemType )

% reads input
xsref = coordsElem( [ 1   7  ] ) ;
ysref = coordsElem( [ 1+2 7+2] ) ;
zsref = coordsElem( [ 1+4 7+4] ) ;

ndofpnode =  6 ;
conecElem = [] ;

% converts global coord displacements to local
if elemType == 2

  xsdef = xsref + dispsElem( [ 1   1+ndofpnode   ] ) ;
  ysdef = ysref + dispsElem( [ 1+2 1+ndofpnode+2 ] ) ;
  zsdef = zsref + dispsElem( [ 1+4 1+ndofpnode+4 ] ) ;
	conecElem = [1 2] ;

  [ ~, ~, ~, rotData ] = elementBeamForces( coordsElem(1:end), ones(5, 1)', [1 1 0], 0, dispsElem, [], [], 0  ) ;
  %~ [ ~, ~, ~, rotData ] = elementBeamForces( coordsElem(1:2:end), ones(7, 1), 0, 0, dispsElem, [], []  ) ;
  
  % global rotation matrix 
  Rr     = rotData{2} ;
  %~ Rr     = eye(3);
  
  locDisp = rotData{1} ;
  ul  = locDisp(1)   ;   tl1 = locDisp(2:4);  tl2 = locDisp(5:7) ;

  titax = [0 0]' ;
  titay = [0 0]' ;
  titaz = [0 0]' ;
  %~ titax = [tl1(1) tl2(1)]' ;
  %~ titay = [tl1(2) tl2(2)]' ;
  %~ titaz = [tl1(3) tl2(3)]' ;

elseif elemType == 3

  xsdefA = coordsElem( [ 1   7  ] ) + dispsElem( [ 1   7  ] ) ;
  ysdefA = coordsElem( [ 1+2 7+2] ) + dispsElem( [ 1+2 7+2] ) ;
  zsdefA = coordsElem( [ 1+4 7+4] ) + dispsElem( [ 1+4 7+4] ) ;

  nPlotPoints    = 10 ; 

  [ ~, ~, ~, rotData ] = elementBeamForces( coordsElem(1:end), ones(5, 1)', [1 1 1], 0, dispsElem, [], [] ,0 ) ;
  
  % 
  locDisp = rotData{1} ;
  ul  = locDisp(1)   ;   tl1 = locDisp(2:4) ;  tl2 = locDisp(5:7) ;

  % global rotation matrix 
  Rr     = rotData{2} ;

  l      = sqrt( sum( (xsref(2)-xsref(1))^2 + (ysref(2)-ysref(1))^2  + (zsref(2)-zsref(1))^2 ) ) ;
  
  xsglo  = linspace( xsdefA(1) , xsdefA(2), nPlotPoints ) ;
  ysglo  = linspace( ysdefA(1) , ysdefA(2), nPlotPoints ) ;
  zsglo  = linspace( zsdefA(1) , zsdefA(2), nPlotPoints ) ;
  
  % loc axial coords to evaluate magnitudes 
  xsloc  = linspace( 0 , l, nPlotPoints )' ;
  
  ux     = zeros( size(xsloc) ) ; 
  uy     = zeros( size(xsloc) ) ; 
  uz     = zeros( size(xsloc) ) ;
  
  titax	 = zeros( size(xsloc) ) ;
  titay	 = zeros( size(xsloc) ) ;
  titaz	 = zeros( size(xsloc) ) ;
  
  LocAxialdofs  = [ 1  1+ndofpnode                        ] ;
  LocTorsidofs  = [ 2  2+ndofpnode                        ] ;
  LocBendXYdofs = [ 3  6           3+ndofpnode 6+ndofpnode] ;
  LocBendXZdofs = [ 5  4           5+ndofpnode 4+ndofpnode] ;

  localUelem = zeros(2,1);
  localUelem( LocAxialdofs  ) = [ 0; ul ] ;
  localUelem( LocTorsidofs  ) = [ tl1(1); tl2(1) ] ;
  localUelem( LocBendXZdofs ) = [ 0; tl1(2); 0; tl2(2) ] ;
  localUelem( LocBendXYdofs ) = [ 0; tl1(3); 0; tl2(3) ] ;

  for i=1:nPlotPoints,
  
    N = bendingInterFuns( xsloc(i), l, 0 );

    Nlin1 = (l - xsloc(i))/l ;
    Nlin2 = xsloc(i)/l       ;

    % bending
    Rb      = eye(4) ;
    Rb(2,2) = -1     ;
    Rb(4,4) = -1     ;

    uy(i)  = N       * localUelem( LocBendXYdofs ) ;
    uz(i)  = N * Rb' * localUelem( LocBendXZdofs ) ;
    
    % stretching
    ux(i) = Nlin1 * localUelem(LocAxialdofs(1)) + Nlin2 * localUelem(LocAxialdofs(2)) ;
		
    % torsion
		titax(i) = Nlin1 * localUelem(LocTorsidofs(1)) + Nlin2 * localUelem(LocTorsidofs(2)) ;
		
    % rots
    Nrots    = bendingInterFuns( xsloc(i), l, 1 ) ;
    titaz(i) = Nrots       * localUelem( LocBendXYdofs ) ; 
    %~ titay(i) = Nrots * Rb' * localUelem( LocBendXZdofs ) ;
    titay(i) = - (Nrots * Rb' * localUelem( LocBendXZdofs ) ) ;
    
		if i < nPlotPoints
			conecElem = [ conecElem ; (i-1)+1 (i-1)+2 ] ;
		end
		
  end
  
  XsLoc = [ xsloc+ux uy uz];
  XsGlo = Rr * XsLoc' ;
  
  xsdef = ( xsdefA(1) + XsGlo(1,:) )' ;
  ysdef = ( ysdefA(1) + XsGlo(2,:) )' ;
  zsdef = ( zsdefA(1) + XsGlo(3,:) )' ;
  
  %~ for k=1:length( titax)
    %~ auxtitaGlo = logar( expon( [ titax(k) titay(k) titaz(k) ]' ) * Rr ) ;
    %~ auxtitaGlo = logar( Rr * expon( [ titax(k) titay(k) titaz(k) ]' ) ) ;
    %~ titax(k) = auxtitaGlo(1);
    %~ titay(k) = auxtitaGlo(2);
    %~ titaz(k) = auxtitaGlo(3);
  %~ end
end
