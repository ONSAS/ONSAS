% ------------------------------------------------------------------------------
% Test problem: beam with springs 
% ------------------------------------------------------------------------------

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

% --------------------------------------------------------------------------------------------------

inputONSASversion = '0.1.8';

problemName = 'beamLSprings' ;

% Constitutive properties

E  = 210e6 ; % Young modulus
nu = 0.3   ; % Poisson
G  = E / ( 2*(1+nu) ) ; % Shear modulus
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu] ;

% Geometrical properties

l1 = 1 ; % length of beam 1
l2 = 1 ;
t = 8e-3 ; % thickness
phi = 100e-3 ; % diameter

% beam 1

A1  = pi*phi^2/4 - pi*(phi-2*t)^2/4 ;     
Iy1 = pi*(phi^4-(phi-2*t)^4)/64 ;
Iz1 = Iy1 ;
It1 = 2*Iz1 ;

% beam 2 

A2  = pi*phi^2/4 - pi*(phi-2*t)^2/4 ;     
Iy2 = pi*(phi^4-(phi-2*t)^4)/64 ;
Iz2 = Iy2 ;
It2 = 2*Iz2 ;


secGeomProps = [ A1 Iy1 Iz1 It1 ; ...
                 A2 Iy2 Iz2 It2 ] ;
                 

% Nodes matrix

Nodes = [ 0 0 0 ; ...
          l2 0 0] ;
          

% Element conectiity

Conec = [ 1 2 0 0 1 1 2 ] ;

% Supports

nodalSprings = [ 1  3*E*Iz1/l1^3  2*E*Iy1/l1^2  E*A1/l1  (G*It1/l1)  3*E*Iy1/l1^3  2*E*Iz1/l1^2 ] ;

% External loads
Fz = -4 ;
nodalConstantLoads   = [ 2  0 0 0 0 Fz 0 ] ;

% Analysis parameters

nonLinearAnalysisBoolean = 0 ; 
printflag = 2;
linearDeformedScaleFactor = 10.0;

% Analytic Solution

analyticSolFlag = 3 ;
G  = E / ( 2*(1+nu) ) ;
Mt = Fz*l2 ;

analytSol = [ Fz*l1^3/(3*E*Iy1)+Fz*l2^3/(3*E*Iy2)+(Mt*l2/(G*It2)) ] ;  

analyticSolDofs = 11 ;
analyticCheckTolerance = 1e-4 ;
% ------------------------------------------------------------------------------
