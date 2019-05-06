% ---------------------------
% Test problem: Cylinder with to internal pressure
% ---------------------------

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

inputONSASversion = '0.1.8'; problemName = 'Cylinder' ;

% Pressure
p = -0.1 ;

% Mesh data
[Nodes, Conec, nodalConstantLoads, nodalSprings] = msh2input( './input/canio.msh', p ) ;

% Material properties
E = 210e9 ; nu = 0.3 ;
  
hyperElasParams = cell(1,1) ;  
hyperElasParams{1} = [ 1 E nu ] ;

% Sections
secGeomProps = [ 0 0 0 0 ] ; Rint = 0.2 ; Rext = 0.24 ;

% Analysis parameters
nonLinearAnalysisBoolean = 0 ; linearDeformedScaleFactor = 10.0 ;

% Plot options
plotParamsVector = [ 3 ] ; printflag = 2 ;

% Analytic solution
analyticSolFlag = 3 ; p = abs(p) ; r = Rext ;
a = ( (1+nu)*(1-2*nu)*Rint^2*p ) / ( E*(Rext^2-Rint^2) ) ;
b = ( (1+nu)*Rint^2*Rext^2*p )   / ( E*(Rext^2-Rint^2) ) ;
analytSol = a*r + b/r ; analyticSolDofs = [ 6*(6-1)+1 ] ;
analyticCheckTolerance = 1e-3 ;
