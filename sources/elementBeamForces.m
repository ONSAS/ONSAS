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

function  [ fs, ks, stress ]= elementBeamForces( ...
  xs, params, booleanCSTangs, solutionMethod, Ue, Udote, Udotdote ) ;

% parameters
E   = params(1) ;
G   = params(2) ;
A   = params(3) ;

Iyy = params(4) ;
Izz = params(5) ;
J   = params(6) ;

rho = params( 7 ) ;

% auxiliar matrices
I3 = eye(3)     ;
O3 = zeros(3)   ;
O1 = zeros(1,3) ;

permutIndxs = [1:2:5 2:2:6 ([1:2:5]+6) ([2:2:6]+6) ] ;

dg       = Ue      ( permutIndxs ) ;
ddotg    = Udote   ( permutIndxs ) ;
ddotdotg = Udotdote( permutIndxs ) ;

% global thetas
tg1 = dg (4:6);
tg2 = dg (10:12);

% rotation matrices
Rg1 = expon( tg1 ) ;
Rg2 = expon( tg2 ) ;

x21 = xs(4:6) - xs(1:3) ;
d21 = dg(7:9) - dg(1:3) ;


lo = sqrt(x21'*x21);
%~ l  = sqrt( (x21+d21)' * (x21+d21) ) ;
%~ l  = norm( x21 + d21 ) ;

l = sqrt( sum( ( x21+d21).^2 ) ) ;

%~ if norm(imag(dg))>0
  %~ u, d21, l, lo, imag(d21), dg
  %~ uimprov = ( l^2 - lo^2 ) / (lo + l)
%~ end

% rotation matrix to reference configuration
Ro = beamRefConfRotMat( x21 ) ;

% --- rigid rotation ---

% deformed x axis
e1 = ( x21 + d21 ) / l ;

q1 = Rg1 * Ro * [0 1 0]' ;
q2 = Rg2 * Ro * [0 1 0]' ;
q  = ( q1 + q2 ) / 2 ;

% deformed z local axis
e3 = cross (e1, q) ;
e3 = e3 / norm(e3) ; % normalization

% deformed z local axis
e2 = cross (e3, e1);

% rotation matrix
Rr = [ e1 e2 e3 ] ;
% -------------------


% --- local displacements ---

% axial strain
u  = l - lo;

% local rotations
Re1 = Rr'*Rg1*Ro;
Re2 = Rr'*Rg2*Ro;

tl1 = logar(Re1);
tl2 = logar(Re2);

locDisp = [ u tl1' tl2' ] ;
% -----------------------


% --- local force vector and tangent stiffness matrix ---
[fl, kl, strain, stress] = beamLocalStaticForces (u, tl1, tl2, lo, E, G, A, Iyy, Izz, J ) ;
% -------------------------------------------------------


q  = Rr' *  q ;
q1 = Rr' * q1 ;

nu = q(1)/q(2);
nu11 = q1(1)/q(2);
nu12 = q1(2)/q(2);
nu21 = 2*nu-nu11;
nu22 = 2-nu12;


% transformation to the new local coordinates

De1 = invTs( tl1 ) ;
De2 = invTs( tl2 ) ;

% matrix for transformation between global and relative rotations/moments
H  = [  1   O1   O1 ; ... 
       O1' De1   O3 ; ...
       O1'  O3  De2 ] ;

fe = H' * fl ;
 %~ [      fl(  1)  
       %~ De1'*fl(2:4)
       %~ De2'*fl(5:7)] ;

Dh1 = dinvTs( tl1, fl(2:4) ) * De1 ;
Dh2 = dinvTs( tl2, fl(5:7) ) * De2 ;

Kh = [ 0   O1   O1
      O1' Dh1   O3
      O1'  O3  Dh2 ] ;

ke = H' * kl * H + Kh ;

% transformation to the global coordinates
r = [-e1;0;0;0;e1;0;0;0];

B = [ r'
   -nu/l*e3' (1-nu12/2)*e1'+nu11/2*e2'  nu/l*e3' 1/2*(-nu22*e1'+nu21*e2')
   -e3'/l e2' e3'/l 0 0 0
    e2'/l e3' -e2'/l 0 0 0
   -nu/l*e3' 1/2*(-nu12*e1'+nu11*e2')  nu/l*e3' (1-nu22/2)*e1'+nu21/2*e2'
   -e3'/l 0 0 0 e3'/l e2'
    e2'/l 0 0 0 -e2'/l e3'];

fg = B' * fe ;

A=(I3-e1*e1')/l;

Dr=[A  O3 -A  O3
    O3 O3  O3 O3
   -A  O3  A  O3
    O3 O3  O3 O3];

G=[0   0    nu/l  nu12/2  -nu11/2  0  0  0    -nu/l  nu22/2  -nu21/2  0
   0   0    1/l     0        0     0  0  0    -1/l     0        0     0
   0  -1/l  0       0        0     0  0  1/l   0       0        0     0]';

II=[O3 I3 O3 O3
    O3 O3 O3 I3];

P=II-[G';G'];

F=P'*fe(2:7);

sF=[skew(F(1:3))
    skew(F(4:6))
    skew(F(7:9))
    skew(F(10:12))];

EE=[Rr O3 O3 O3
    O3 Rr O3 O3
    O3 O3 Rr O3
    O3 O3 O3 Rr];

nab=[0
    (nu*(fe(2)+fe(5))+fe(3)+fe(6))/l
    (fe(4)+fe(7))/l]; 

Kg = B' * ke * B + Dr * fe(1) - EE*sF*G'*EE' + EE*G*nab*r' ;


% --- transformation to the new global coordinates ---

Dg1=Ts(tg1);
Dg2=Ts(tg2);

q=[fg(1:3)
   Dg1'*fg(4:6)
   fg(7:9)
   Dg2'*fg(10:12)];

Dk1=dTs(tg1,fg(4:6));
Dk2=dTs(tg2,fg(10:12));

H=[I3 O3  O3 O3
   O3 Dg1 O3 O3
   O3 O3  I3 O3
   O3 O3  O3 Dg2];

Kt = H' * Kg * H ;

Kt( 4:6 , 4:6 ) = Kt( 4:6 , 4:6 ) + Dk1 ;
Kt(10:12,10:12) = Kt(10:12,10:12) + Dk2 ;

%~ Kt = (Kt+Kt')/2;

Finte   = zeros(size(q)) ;
dofscomb = [ 1:2:5 2:2:6 7:2:11 8:2:12 ] ;

Finte( dofscomb ) = q ;
KTe = zeros( size(Kt));

if booleanCSTangs == 1

  step = 1e-4 * norm(x) ;
  
  for i=1:12
    ei = zeros(12,1);   ei(i) = j ;
    
    FinteComp = elementBeamInternLoads( x, dg + ei*step, params, 0 ) ;
    
    KTe(:,i) = imag( FinteComp ) / step;
    
    %~ if i==1
      %~ holaaafintecomp = FinteComp(1) ;    
    %~ ei
%~ FinteComp
%~ stop
    %~ end
  end
  %~ KTeCS = KTe ;
%~ KTe = zeros( size(Kt));
  %~ KTe( dofscomb, dofscomb ) = Kt ;
  %~ normareldif = norm( KTeCS - KTe ) / norm( KTe )
  %~ dife = KTeCS - KTe
  %~ normareldif11 = norm( KTeCS(1,1) - KTe(1,1) ) / norm( KTe(1,1) )
  %~ entridif = [ KTeCS(1,1) KTe(1,1) holaaafintecomp ]
  %~ holacomplejos = [ KTeCS(1,1) holaaafintecomp ]

  %~ full(dife)
    
  %~ stop
else

  KTe( dofscomb, dofscomb ) = Kt ;
end

fs = {Finte} ;
ks = {KTe};





