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
function  [ fs, ks, stress, rotData ] = frame_internal_force( ...
  elemCoords, elemCrossSecParams, elemConstitutiveParams, Ue ) ;
  
  % element coordinates
  xs = elemCoords(:) ;

  % ----- material and geometric params ------
  E   = elemConstitutiveParams(2) ;
  nu  = elemConstitutiveParams(3) ;
  G   = E/(2*(1+nu))              ;
  % ----- extract cross section properties ---
  [Area, J, Iyy, Izz, ~] = crossSectionProps ( elemCrossSecParams, 0 ) ; % select a ficticious elemrho 
  % ------------------------------------------

  % compute corotational matrices rotation 
  [R0, Rr, Rg1, Rg2, Rroof1, Rroof2] = corotRotMatrices( Ue, elemCoords ) ;

  % --- global displacements ---
  % permut indexes according to Battini's nomenclature
  dg = switchToBattiniNom( Ue ) ;
  % global thetas
  tg1 = dg(  4:6  ) ;
  tg2 = dg( 10:12 ) ;
  % -------------------------------

  % length and coords of the element  
  [x21, d21, l, l0] = corotLenCoords(xs ,dg) ;

  % --- local displacements ---

  % axial displacement
  u   = l - l0 ;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % temporary
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  global pretension_strain
  if ~isempty( pretension_strain )
    % axial displacement
    u   = u + pretension_strain * l0;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % local rotations
  tl1 = logar( Rroof1 ) ;
  tl2 = logar( Rroof2 ) ;  
  locDisp = [ u tl1' tl2' ] ;
  % -------------------------------

  % --- auxiliary vector and matrices  ---
  % aux zero matrices 
  [I3, O3, O1, II] = corotZeros() ;

  % auxiliary q base 
  q1g = Rg1 * R0 * [0 1 0]' ;
  q2g = Rg2 * R0 * [0 1 0]' ;
  qg  = ( q1g + q2g ) / 2   ;

  [nu, nu11, nu12, nu21, nu22, e1, e2, e3, r, Gaux, P, EE ] = corotVecMatAuxStatic(...
                                                                R0, Rr, Rg1, Rg2, l, II, O3, O1);
  % -------------------------------

  % --- local force vector and tangent stiffness matrix ---
  [fl, kl, strain, stress] = beamLocalStaticForces (u, tl1, tl2, l0, E, G, Area, Iyy, Izz, J ) ;

  % -- term to equal virtual work of moments in local an global --
  % transformation to the new local coordinates
  De1 = invTs( tl1 ) ;
  De2 = invTs( tl2 ) ;

  % matrix for transformation between global and relative rotations/moments
  H  = [  1   O1   O1 ; ...
         O1' De1   O3 ; ...
         O1'  O3  De2 ] ;

  fe = H' * fl ;

  Dh1 = dinvTs( tl1, fl(2:4) ) * De1 ;
  Dh2 = dinvTs( tl2, fl(5:7) ) * De2 ;

  Kh = [ 0   O1   O1
        O1' Dh1   O3
        O1'  O3  Dh2 ] ;

  ke = H' * kl * H + Kh ;
  % -------------------------------------------------------------

  %  -------------transformation to the global coordinates-------
  B = [ r'
     -nu/l*e3' (1-nu12/2)*e1'+nu11/2*e2'  nu/l*e3' 1/2*(-nu22*e1'+nu21*e2')
     -e3'/l e2' e3'/l 0 0 0
      e2'/l e3' -e2'/l 0 0 0
     -nu/l*e3' 1/2*(-nu12*e1'+nu11*e2')  nu/l*e3' (1-nu22/2)*e1'+nu21/2*e2'
     -e3'/l 0 0 0 e3'/l e2'
      e2'/l 0 0 0 -e2'/l e3' ];

  fg = B' * fe ;

  A  = (I3-e1*e1')/l;

  Dr = [A  O3 -A  O3
        O3 O3  O3 O3
       -A  O3  A  O3
        O3 O3  O3 O3];

  F = P' * fe(2:7);

  sF=[ skew( F( 1:3 )  )
       skew( F( 4:6 )  )
       skew( F( 7:9 )  )
       skew( F( 10:12) ) ] ;


  nab=[ 0
        ( nu * ( fe( 2 ) + fe( 5 ) ) + fe( 3 ) + fe( 6 ) ) / l
        ( fe( 4 ) + fe( 7 ) ) / l ] ;

  Kg = B' * ke * B + Dr * fe(1) - EE * sF * Gaux' * EE' + EE * Gaux * nab * r' ;

  Dg1 = Ts( tg1 ) ;
  Dg2 = Ts( tg2 ) ;

  q=[ fg( 1:3 )
      Dg1' * fg( 4:6 )
      fg( 7:9 )
      Dg2' * fg( 10:12 ) ] ;

  Dk1 = dTs( tg1, fg( 4:6 ) )   ;
  Dk2 = dTs( tg2, fg( 10:12 ) ) ;

  H=[ I3 O3  O3 O3 
      O3 Dg1 O3 O3 
      O3 O3  I3 O3 
      O3 O3  O3 Dg2 ];

  Kt = H' * Kg * H ;

  Kt(  4:6 , 4:6  ) = Kt(  4:6 , 4:6  ) + Dk1 ;
  Kt( 10:12,10:12 ) = Kt( 10:12,10:12 ) + Dk2 ;

  % make the internal tangent matrix symmetric to improve convergence order
  Kt = ( Kt + Kt' ) / 2;

  % --- Write it back to ONSAS nomenclature [fx mx,....]  ---
  Finte = switchToONSASNom(q) ;

  KTe = zeros( size(Kt) );


  % boolean to compute the tangents with complex-step
  booleanCSTangs = 0 ;
  if booleanCSTangs == 1 

    step = 1e-4 * norm(x) ;

    for i=1:12
      ei = zeros(12,1);   ei(i) = j ;

      FinteComp = elementBeamInternLoads( x, dg + ei*step, params, 0 ) ;

      KTe(:,i) = imag( FinteComp ) / step;
    end

  else

    KTe = switchToONSASNom(Kt) ;

  end
  % -------------------------------------------------------

  % -------- Compress output -----------
  fs = {Finte} ;
  ks = {KTe};
  rotData = {locDisp, Rr} ;
  % ------------------------------------