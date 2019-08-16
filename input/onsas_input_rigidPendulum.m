% TEST example springmass
% ------------------------------------

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


% auxiliar numerical data
Es = 100 ;
A  = 1 ;
l0 = 2 ;

rhoprob =  1    ;
ceda    =  0.02 ;
nu      =  0    ;

% method
timeIncr   =  0.01    ;
finalTime  =  2       ;
nLoadSteps = finalTime/timeIncr ;

DeltaNW    =  0.5               ;
AlphaNW    =  0.25              ;

% tolerances
stopTolDeltau = 1e-10           ; 
stopTolForces = 1e-10           ;
stopTolIts    = 30            ;
% ------------------------------------


% --- general data ---
inputONSASversion = '0.1.8';

problemName = 'rigidPendulum' ;
% ------------------------------------

% --- structural properties ---
rho = rhoprob ;
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 Es nu rho] ;


secGeomProps = [ A 0 0 0 ] ;

nodalSprings = [ 1  inf  0  inf  0  inf 0   ...
               ];

Nodes = [    0  0  0 ; ...
            l0  0  0 ] ;


Conec = [ 1 2 0 0 1 1 1 ] ; 

%~ loadFactorsFunc = @(t) p0 *sin( omegaBar*t) ; 

% -------------------
%~ nodalVariableLoads   = [ 2  1  0  0  0  0  0 ];
% or
%~ nodalVariableLoads   = [ 2  0  0  0  0  0  0 ];
nodalConstantLoads   = [ 2  0  0  0  0  -rhoprob*l0*A*.5*9.81  0 ];
% -------------------

controlDofInfo = [ 2 5 -1 ] ;
% ------------------------------

% ------------------------------
% analysis parameters
dynamicAnalysisBoolean   = 1 ; 

% initial conditions
u0   = 0;
udot0= 0;

%~ nonHomogeneousInitialCondU0 = [ 2 1 u0 ] ;

%~ nonHomogeneousInitialCondUdot0 =[ 2 1 udot0 ] ;

numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

plotParamsVector = [2 5 ];
printflag = 2 ;

analyticSolFlag = 0 ;
% ------------------------------------------------
