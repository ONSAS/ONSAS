% ------------------------------------------------------------------------------
% Test problem: Plate submitted to distributed load on surface with beams on the edges and columns 
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
problemName = 'PlateBeamColumn' ;

ndofpnode = 6 ;

% Material properties
E = 8500*(25+8)^(1/3)*100 ;
nu  = 0 ;

hyperElasParams = cell(1,1) ;  
hyperElasParams{1} = [ 1 E nu ] ;

% Geometric properties

% Slab
t = 0.09 ; 
L = [ 4 4 ] ;

% Rigid beam
a1 = 0.9 ; 
b1 = 0.2 ; 
% Flexible beam
a2 = 0.5 ;
b2 = 0.5 ;
% Columns
a3 = 0.3 ;
b3 = 0.13 ; 

constRigid = 1e4 ;
constColumns = 1e4 ;
constFlexible = 1e-9 ;

%  rigid beams
Ab = a1*b1 ;
Izb = a1*b1^3/12 ;
Iyb = b1*a1^3/12*constRigid ;
beta = 0.141 ;
Itb = beta * b1 * a1^3*constRigid ;

% flexible beams
Ab2 = a3*b3 ;
Izb2 = a3*b3^3/12 ;
Iyb2 = b3*a3^3/12*constFlexible ;
beta = 0.141 ;
Itb2 = beta * b3 * a3^3*constFlexible ;

% columns
Ac = a2*b2*constColumns ;
Izc = a2*b2^3/12*constColumns ;
Iyc = b2*a2^3/12*constColumns ;
beta = 0.141 ;
Itc = beta * b2 * a2^3*constColumns ;

secGeomProps = [ Ab   Iyb   Izb   Itb   ; ...
                 Ab2  Iyb2  Izb2  Itb2  ; ...
                 Ac   Iyc   Izc   Itc   ] ;

% Mesh
nx  = 18 ;                    % # divitions in direction "x"
ny      = ceil(L(2)/L(1)*nx) ;  % # divitions in direction "y"
nel     = [ nx ny ]' ;
nnos    = nel + 1 ;  % # of nodes in each direction

% Nodes
lins1   = linspace( 0 , L(1) , nnos(1) )' ;
lins2   = linspace( 0 , L(2) , nnos(2) )' ;

Nodes = [ ] ;  
for i = 1:nnos(2) 
  Nodes( ( (nnos(1) * (i-1) + 1) : nnos(1)*i ) ,:) = [ lins1 lins2(i)*ones(nnos(1),1) zeros(nnos(1),1) ] ;
end

nnodesplate = size(Nodes,1) ;
midnodeplate = floor(nnodesplate/2) +1 ;
h = 4 ;

Nodes = [ Nodes ; ...
          0 0 -h ; ... 
          4 0 -h ; ...
          4 4 -h ; ...
          0 4 -h ] ;

nnodes = size(Nodes,1) ;

% Conectivity matrix 
Conec = [ ] ; 

% Plates
for j = 1:nel(2)
  for i = 1:nel(1)
    intri = (i-1)+1+(j-1)*nel(1) ;
    Conec( intri , : ) = [(j-1)*nnos(1)+i    (j-1)*nnos(1)+i+1   j*nnos(1)+i+1   j*nnos(1)+i] ;
  end
end
Conec = [ Conec ones(size(Conec,1),1) zeros(size(Conec,1),1) ones(size(Conec,1),1)*4 ] ;
% Beams
auxab = [] ;
for i = 1:nx
  auxab = [ auxab ; i i+1 0 0 1 2 2 ] ;
end
auxar = [] ;
for i = 1:nx
  auxar = [ auxar ; (nx+1)*ny+i (nx+1)*ny+i+1 0 0 1 2 2 ] ;
end
auxizq = [] ;
for i = 1:ny
  auxizq = [ auxizq ; (i-1)*nnos(1)+1 (i)*nnos(1)+1 0 0 1 1 2 ] ;
end
auxder = [] ; 
for i = 1:ny
  auxder = [ auxder ; i*nnos(1) (i+1)*nnos(1) 0 0 1 1 2 ] ;
end
% Columns
auxcolumn = [ 1 nnodesplate+1            0 0 1 3 2 ; ...
              nnos(1) nnodesplate+2      0 0 1 3 2 ; ...
              nnos(1)*ny+1 nnodes        0 0 1 3 2 ; ...
              nnodesplate nnodesplate+3  0 0 1 3 2 ] ; 
               
Conec = [ Conec ; auxab ; auxder ; auxar ; auxizq ; auxcolumn ] ; 

nelems = size(Conec,1) ; 

% Boundary conditions
nodalSprings = [ nnodesplate+1 inf 0 inf 0 inf 0 ; ...
                 nnodesplate+2 inf 0 inf 0 inf 0 ; ...
                 nnodesplate+3 inf 0 inf 0 inf 0 ; ...
                 nnodesplate+4 inf 0 inf 0 inf 0 ] ;

% Nodal Loads
q = -0.15 ; % ton/m2
Lx = L(1) ;
Ly = L(2) ;
a = L(1)/(nx*2) ;
b = L(2)/(ny*2) ;

elemsPlate = unique(find(Conec(:,7)==4)) ;

unifLoad = [ (elemsPlate) ones(length(elemsPlate),1) zeros(length(elemsPlate),2) ones(length(elemsPlate),1)*q ] ;

% Analysis parameters
nonLinearAnalysisBoolean = 0 ; 
printflag = 2 ;
linearDeformedScaleFactor = 10.0 ;

% Plot options
plotParamsVector = [ 3 ] ;
printflag = 2 ;

% Analytic solution 

analyticSolFlag = 3 ;

D = E * t^3 / ( 12 * (1-nu^2) ) ;
dof = midnodeplate*ndofpnode-1 ;
analytSol = q * (L(1)/2)^2*(L(1)-L(1)/2)^2/(24*D) ;
analyticSolDofs = [dof] ;
analyticCheckTolerance = 1e-3 ;
