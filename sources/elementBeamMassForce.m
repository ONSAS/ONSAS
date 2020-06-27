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


function [Fine,MassMatrix,GyroMatrix] = elementBeamMassForce(xelem, Dte, Ddote, Ddotdote, params,Jrho )

% ----- Material and geometrical params -----
E   = params(1) ;
G   = params(2) ;
A   = params(3) ;

Iyy = params(4);
Izz = params(5);
J   = params(6);

rho = params(7) ;

%Irho es el tensor de inercia de rotacional de segundo orden en la config indeformada

%---------Global displacements and volocities---------
% permutation
permutIndxs = [1:2:5 2:2:6 ([1:2:5]+6) ([2:2:6]+6) ] ;

dg       = Dte     ( permutIndxs ) ;
ddotg    = Ddote   ( permutIndxs ) ;
ddotdotg = Ddotdote( permutIndxs ) ;
% ----------------------------------------------------

%---------Rotations matrix Rg Ro Rr Rtecho---------

%Global rotation Rg
tg1 = dg(  4:6  ) ;
tg2 = dg( 10:12 ) ;% tita global nodo 1 y 2 

Rg1 = expon(tg1); Rg2 = expon(tg2); %Matrices Rg de rotacion globales

x21 = xelem(4:6) - xelem(1:3); %vector  que va del nodo 2 a 1 (en condiciones indeformadas)

I3 = eye(3);

d21 = dg(7:9) - dg(1:3); %Vector desplzamiento que va de 2 a 1

lo = sqrt( ( x21       )' * ( x21       ) ) ; %longitud inicial
l  = sqrt( ( x21 + d21 )' * ( x21 + d21 ) ) ; % largo deformado 
u  = l-lo; %variable del estiramiento axial

%Material rotation Ro
Ro = rotRol( x21 ) ; %Matriz que va del material a la comfgiguracion indeformada

%Rigid rotation Rr

r1 = (x21+d21)/l;%vector que une el punto 2 con el 1 en la deformada
q1 = Rg1*Ro*[0;1;0];%vector q1 y q2
q2 = Rg2*Ro*[0;1;0];
q  = (q1+q2)/2; %vector q muestra en cierto sentido el angulo de torsión

r3 = cross (r1, q);

r3 = r3 / norm( r3 ) ;% calcula e3 como en el paper

r2 = cross (r3, r1);%e2 y ya queda definida la matriz

Rr = [r1 r2 r3];

q    = Rr'*q; % Parametros para la construccion de G
q1   = Rr'*q1;
nu   = q(1)/q(2);
nu11 = q1(1)/q(2);
nu12 = q1(2)/q(2);
nu21 = 2*nu-nu11;
nu22 = 2-nu12;


%cuelga la transofrmada de los vectores y eso resulta la matriz que da el movimiento deformable
%r1 es e1 ene el paper
%Rotation Rroof

Re1 = Rr'*Rg1*Ro;%Ec 18
Re2 = Rr'*Rg2*Ro;
tl1 = logar(Re1);%Ec 19
tl2 = logar(Re2);

locDisp = [ u; tl1; tl2 ];

%---------Local interplation ---------
N1 = @(x) 1 -x/lo			 ;
N2 = @(x) x/lo				 ;
N3 = @(x) x*(1-x/lo)^2		 ;
N4 = @(x) -(1-x/lo)*(x^2)/lo ; 
N5 = @(x) (1-3*x/lo)*(1-x/lo);
N6 = @(x) (3*x/lo-2)*(x/lo)	 ;

N7 = @(x) N3(x)+N4(x)		 ;
N8 = @(x) N5(x)+N6(x)-1		 ;

%---------Compute kinematic variables ---------

%Matrix to compute udot y udtotdot
G  = [0   0    nu/l  nu12/2  -nu11/2  0  0  0    -nu/l  nu22/2  -nu21/2  0
      0   0    1/l     0        0     0  0  0    -1/l     0        0     0
      0  -1/l  0       0        0     0  0  1/l   0       0        0     0]';

O3 = zeros(3,3);
O1 = zeros(1,3);
I3 = eye(3);


II = [O3 I3 O3 O3
	  O3 O3 O3 I3]; 

P  =  II-[G';G']; %Ec 32

P1    = @(x) [0   0       0     0   0    0   ;
			  0   0     N3(x)   0   0   N4(x);
			  0 -N3(x)	  0    0 -N4(x)  0	]; %Ec	38

P2    = @(x) [ N1(x)  0     0    N2(x)  0      0;
				0   N5(x)   0     0     N6(x)  0;
				0    0     N5(x)  0      0    N6(x)]; %Ec 39

EE     = [Rr O3 O3 O3
		  O3 Rr O3 O3
		  O3 O3 Rr O3
		  O3 O3 O3 Rr]; %Ec 30 es E en el texto
     
ul    = @(x)P1(x)*[tl1 ;
					tl2 ];%Ec 38

N     = @(x) [N1(x)*I3 O3 N2(x)*I3 O3]; %falta las dimensiones 

H1    = @(x) N(x) +P1(x)*P-skew(ul(x))*G';% Ec 59

wdoter= G'*EE'*ddotg;%Ec 65

r  	  = [-r1' O1 r1' O1]; % Ec A.5
A1    = [   O1      O1     O1       O1;
		  0 -1  0   O1 	 0  1  0	O1;
		  0  0 	-1  O1	 0  0  1    O1]; %Ec A.4
	  
udotl = @(x)  P1(x)*P*EE'*ddotg; %Ec A.9

H1dot = @(x)  N7(x)/(l^2)*A1*(r*ddotg)-skew(udotl(x))*G'; %Ec A.8

ET = [skew(wdoter)      O3         		O3 			 O3			;
			O3		skew(wdoter)  		O3   		 O3			;		
			O3			O3 			skew(wdoter)     O3			;							
			O3			O3  			O3       skew(wdoter)   ];

C1 = @(x)  skew(wdoter)*H1(x) + H1dot(x) -H1(x)*ET; % Ec  66

udot    = @(x) Rr*H1(x)*EE'*ddotg; %Ec 61
udotdot = @(x) Rr*H1(x)*EE'*ddotdotg+Rr*C1(x)*EE'*ddotg; % Ec 67



%Matrix to compute wdot y wdtotdot

H2 = @(x) P2(x)*P+G'; %Ec 72 se puede usar para comprobar con ec A.10

wdot  = @(x) Rr*H2(x)*EE'*ddotg;%Ec74


A2    = [   O1    O1    O1     O1;
		  0 0  1  O1  0 0 -1   O1;
		  0 -1 0  O1  0 1  0   O1];%Ec A.12
		  
H2dot    = @(x)	N8(x)/l^2*A2*(r*ddotg) ;%Ec A.14	  

C2       = @(x) skew(wdoter)*H2(x) + H2dot(x) - H2(x)*ET ;%Ec 76	  

wdotdot  = @(x) Rr*H2(x)*EE'*ddotdotg  + Rr*C2(x)*EE'*ddotg ;%Ec 77

%---------Tensor dyadc of Intertia ---------
%compute Rg(x)
thethaRoof  = @(x) P2(x)*[tl1;tl2];% Ec 39
Rex         = @(x) expon(thethaRoof(x)); %Ec 19 elevado en ambos lados
Rgx  	    = @(x) Rr*Rex(x)*Ro'; 


Irho		= @(x) Rgx(x)*Ro*(Jrho)*(Rgx(x)*Ro)'; %Ec 45
Irhoe       = @(x) Rr'*Irho(x)*Rr;   	 		%Ec 80

%---------Compute interial force by quadrature ---------

InterPoints  = [-sqrt(3/5)   0    sqrt(3/5)];
WhightPoints = [    5/9	   8/9     5/9    ];

IntegrandoForce  = @(x) H1(x)'*Rr'*A*rho*udotdot(x)+H2(x)'*Rr'*(Irho(x)*wdotdot(x)...
									+skew(wdot(x))*Irho(x)*wdot(x));  %Ec 78
%~ irho=Irho(sqrt(3/5))
%~ termino=H2(1)'*Rr'*(Irho(1)*wdotdot(1)+skew(wdot(1))*Irho(1)*wdot(1))
sumForce   = zeros(12,1);
for index = 1:size(InterPoints,2)
	
	sumForce = sumForce +  WhightPoints(index)*IntegrandoForce(l/2*InterPoints(index)+l/2);
   if index == size(InterPoints,2)
		sumForce = sumForce*l/2;
   end	
end
Fine = EE*sumForce;
%---------------- Mass matrix  -------------------


sumMass = zeros(12,12) ;

IntegrandoMassMatrix  = @(x) H1(x)'*A*rho*H1(x)+H2(x)'*Irhoe(x)*H2(x); %Ec87

for index = 1:size(InterPoints,2)
	sumMass = sumMass +  WhightPoints(index)*IntegrandoMassMatrix(l/2*InterPoints(index)+l/2);
   if index == size(InterPoints,2)
		sumMass = sumMass*l/2; %sale de la integracion con  cuadrature
   end	
end

MassMatrix = EE*sumMass*EE';

%---------------- Gyroscopic matrix  -------------------

%~ %Compute C3 and C4

h1 =@(x) H1(x)*ddotg; %Ec B6
h2 =@(x) H2(x)*ddotg;

rElem = [ [-1 0 0]   O1  [1 0 0] O1]; %Ec B10

F1    = [skew(udot(0))' skew(wdot(0))' skew(udot(lo))' skew(wdot(lo))']'; %Chequear con los nodales
%~ F1aux    = [skew(ddotg(1:3))' skew(ddotg(4:6))' skew(ddotg(7:9))' skew(ddotg(10:12))']' %Chequear con los nodales

C3  = @(x) -skew(h1(x))*G'  + (N7(x)/l^2)*A1*(ddotg*rElem)...	
							+skew(wdoter)*P1(x)*P + H1(x)*F1*G'; % B13

C4  = @(x) -skew(h2(x))*G' + (N8(x)/l^2)*A2*ddotg*rElem + H2(x)*F1*G'; %B14

%Compute Gyroscopic Matrix
%~ Irhoe(l)
%~ c1prueba = C1(l/2)
%~ c3prueba =C3(l/2)
IntegrandoGyroMatrix  = @(x) H2(x)'*( ( skew(wdoter)*Irhoe(x) ) - Irhoe(x) * skew(wdoter) ) * H2(x) ...
									+	H1(x)'*A*rho*(C1(x) + C3(x))  + H2(x)'*Irhoe(x)*(C2(x)+C4(x)) ; %Ec88
sumGyro = zeros (12,12);
  for index = 1:size(InterPoints,2)
	sumGyro = sumGyro +  WhightPoints(index)*IntegrandoGyroMatrix(l/2*InterPoints(index)+l/2);
	if index == size(InterPoints,2)
	sumGyro = sumGyro*l/2;
  end

GyroMatrix = EE*sumGyro*EE';


%Cambio de Bases

MassMatrix = Cambio_Base(MassMatrix); % En formato [u1 theta1 u2 theta2 u3 theta3];
GyroMatrix = Cambio_Base(GyroMatrix); % En formato [u1 theta1 u2 theta2 u3 theta3];
Fine       = Cambio_Base(Fine); % En formato [f1 m1 ...];

%Usar la funcion cambio de base 
end


%~ function quadSum = integr( hola )

