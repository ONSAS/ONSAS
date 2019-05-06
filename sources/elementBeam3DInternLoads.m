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

% function for computation of nodal forces and tangent stiffness matrix for 3D co-rotational beam element. Based on files provided by Prof. Jean-Marc Battini.

function [ Finte, KTe, strain, stress, locDisp, Rr] = elementBeam3DInternLoads( x, Ue, params )

E   = params(1);
G   = params(2);
A   = params(3);
Iyy = params(4);
Izz = params(5);
J   = params(6);

p = zeros(size(Ue));
p([1:3]  ) = Ue([1:2:5]  ) ;
p([1:3]+6) = Ue([1:2:5]+6) ;
p([4:6]  ) = Ue([2:2:6]  ) ;
p([4:6]+6) = Ue([2:2:6]+6) ;


% tita global nodo 1
tg1=p(4:6);  tg2=p(10:12);

Rg1 = expon(tg1); Rg2 = expon(tg2);

x21 = x(4:6) - x(1:3);

I3 = eye(3);

d21 = p(7:9) - p(1:3);

lo = sqrt(x21'*x21);
l  = sqrt((x21+d21)'*(x21+d21));
u  = l-lo;

Ro = rotRo1(x21);

% rigid rotation

e1 = (x21+d21)/l;
q1=Rg1*Ro*[0;1;0];
q2=Rg2*Ro*[0;1;0];
q=(q1+q2)/2;

e3= cross (e1, q);

e3=e3/norm(e3);

e2 = cross (e3, e1);

Rr=[e1 e2 e3];

q  = Rr'*q;
q1 = Rr'*q1;
nu = q(1)/q(2);
nu11 = q1(1)/q(2);
nu12 = q1(2)/q(2);
nu21 = 2*nu-nu11;
nu22 = 2-nu12;


%local rotations

Re1 = Rr'*Rg1*Ro;
Re2 = Rr'*Rg2*Ro;
tl1 = logar(Re1);
tl2 = logar(Re2);

locDisp = [ u; tl1; tl2 ];

%local force vector and tangent stiffness matrix

[fl,kl, strain, stress] = b3dJor (u,tl1,tl2,lo, E,G,A,Iyy,Izz,J);


% transformation to the new local coordinates

De1 = invTs(tl1);
De2 = invTs(tl2);

fe=[fl(1)
    De1'*fl(2:4)
    De2'*fl(5:7)];

O3=zeros(3,3);
O1=zeros(1,3);

H=[1   O1  O1 
   O1' De1 O3
   O1' O3  De2];

Dh1 = dinvTs(tl1,fl(2:4))*De1;
Dh2 = dinvTs(tl2,fl(5:7))*De2;

Kh = [0   O1  O1
      O1' Dh1 O3
      O1' O3  Dh2 ] ;

ke = H' * kl * H + Kh;

% transformation to the global coordinates

r=[-e1;0;0;0;e1;0;0;0];

B=[ r'
   -nu/l*e3' (1-nu12/2)*e1'+nu11/2*e2'  nu/l*e3' 1/2*(-nu22*e1'+nu21*e2')
   -e3'/l e2' e3'/l 0 0 0
    e2'/l e3' -e2'/l 0 0 0
   -nu/l*e3' 1/2*(-nu12*e1'+nu11*e2')  nu/l*e3' (1-nu22/2)*e1'+nu21/2*e2'
   -e3'/l 0 0 0 e3'/l e2'
    e2'/l 0 0 0 -e2'/l e3'];

fg=B'*fe;

I3=eye(3);

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
KTe( dofscomb, dofscomb ) = Kt ;

% ==============================================================================
% ==============================================================================
function [R]=expon(t);

al=norm(t);
I=eye(3,3);

if al==0
  R=I;
else
  Rsk=skew(t);
  R=I+sin(al)/al*Rsk+2*(sin(al/2)/al)^2*Rsk^2;
end


% ==============================================================================
% ==============================================================================
function Ro =rotRo1(x);

r1 = x / norm(x) ;

q = [0;0;1];

r2 = cross (q,r1);

if norm(r2)~=0
  r2=r2/norm(r2);
  [r3]=cross(r1,r2);
else
  q=[0;1;0];
  [r3]=cross(r1,q);
  r3=r3/norm(r3);
  [r2]=cross(r3,r1);
end
Ro=[r1 r2 r3];




% ==============================================================================
% ==============================================================================
function [fl,kl, strain, stress] = b3dJor (u, tl1, tl2, L, E, G, A, Iyy, Izz, J);

ltx1 = tl1(1);
lty1 = tl1(2);
ltz1 = tl1(3);

ltx2 = tl2(1);
lty2 = tl2(2);
ltz2 = tl2(3);

% --- internal forces vector --

fl = zeros(7,1);

% Auxiliar epsilon x
strain = u/L ;
% Auxiliar sigma x
stress = E*strain ;

fl(1) = E*A*u/L;

fl(2) = - G * J / L * (ltx2-ltx1) ;
fl(5) = -fl(2);

fl(3) = 2*E*Iyy/L * ( 2*lty1 +   lty2 ) ;
fl(6) = 2*E*Iyy/L * (   lty1 + 2*lty2 ) ;

fl(4) = 2*E*Izz/L * ( 2*ltz1 +   ltz2 ) ;
fl(7) = 2*E*Izz/L * (   ltz1 + 2*ltz2 ) ;

% stiffness matrix
kl = zeros(7,7);

kl(1,1) = E*A/L ;

kl([2 5], [2 5] ) = G * J / L  * [ 1 -1 ; -1 1 ];

kl([3 6], [3 6] ) = 2* E * Iyy / L * [ 2  1 ; 1 2 ];

kl([4 7], [4 7] ) = 2* E * Izz / L * [ 2  1 ; 1 2 ];


% ==============================================================================
% ==============================================================================
function [t]=logar(R);

u=[R(3,2)-R(2,3)
   R(1,3)-R(3,1)
   R(2,1)-R(1,2)];

nu=norm(u);

if nu==0
  t=[0;0;0];
else
  t=asin(nu/2)/nu*u;
end


% ==============================================================================
% ==============================================================================
function [De]=invTs(t);

nt=norm(t);
I=eye(3,3);

if nt==0
  De=I;
else
  b=nt/2;
  a=(sin(b)-b*cos(b))/(nt^2*sin(b));
  M=skew(t);
  De=I-1/2*M+a*M*M;
end





% ==============================================================================
% ==============================================================================
function [Dk]=dTs(t,v);

nt=norm(t);

if nt==0
  Dk=1/2*skew(v);
else
  e=t/nt;
  ev=cross(e,v);

  a1=(cos(nt)-sin(nt)/nt)/nt;
  M1=v*e'-(e'*v)*e*e';

  a2=(1-sin(nt)/nt)/nt;
  M2=e*v'-2*(e'*v)*e*e'+(e'*v)*eye(3);

  a3=sin(nt)/nt-(2*sin(nt/2)/nt)^2;
  M3=ev*e';

  a4=2*(sin(nt/2)/nt)^2;
  M4=skew(v);

  Dk=a1*M1+a2*M2-a3*M3+a4*M4;
end


% ==============================================================================
% ==============================================================================
function [Dh]=dinvTs(t,v);

nt=norm(t);

if nt==0
  Dh=-1/2*skew(v);
else
  a=nt/2;
  eta=(sin(a)-a*cos(a))/(nt^2*sin(a));
  miu=(nt*(nt+sin(nt))-8*sin(a)^2)/(4*nt^4*sin(a)^2);
  I3=eye(3);
  M=skew(t);
  M1=skew(v);
  M2=t*v'-2*v*t'+(t'*v)*I3;
  M3=M*M*v*t';
  Dh=eta*M2+miu*M3-1/2*M1;
end

% ==============================================================================
% ==============================================================================
function [sk] = skew(x);

sk = [   0 -x(3)  x(2) ;
      x(3)    0  -x(1) ;
     -x(2)  x(1)     0 ];
  

% ==============================================================================
% ==============================================================================
function [Dg]=Ts(t);

nt=norm(t);
I=eye(3,3);

if nt==0
  Dg=I;
else
  a=2*(sin(nt/2)/nt)^2;
  b=(1-sin(nt)/nt)/nt^2;
  M=skew(t);
  Dg=I+a*M+b*M*M;
end
