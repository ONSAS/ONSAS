% Copyright (C) 2020, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
%   Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini, Sebastian Toro  
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


function solid3DlinearAnalysis

close all

% =========  1 - Parametros =============================================

% elastic parameters
E        = 4.9e9 ;
nu       = 0.4   ;

% calcula params  
mu    = E / (2.0 * ( 1.0 + nu ) ) ;
Bulk  = E / (3.0 * ( 1.0 - 2.0 * nu ) ) ;

eyetres      = eye(3) ;
eyevoig      = zeros(6,1) ;
eyevoig(1:3) = 1.0 ;


% tensor constitutivo 
ConsMat = zeros(6,6) ;

ConsMat (1:3,1:3) = Bulk  + 2* mu * ( eyetres - 1.0/3.0 ) ;
ConsMat (4:6,4:6) =            mu *   eyetres ;

vectorCarga = [ 1 0 0 ]';
% ======================================================================



% =========  2 - Mallado =============================================
[nodes, Cellsconec, nodoscargados, nodosFijos, fext] = mallayahecha ;

nodoscargados
nodosFijos

nnod     = size(nodes,1) ;
nodesdef = nodes ;
ntet     = size(Cellsconec,1) 

DiriDofs      = nods2dofs( nodosFijos , 3 );
% ======================================================================================

DiriDofs


% = 3 - Preproceso ======================================================================
NeumDofs = (1: 3*nnod ) ;
NeumDofs( DiriDofs ) = [] ;

DIRECTA = 1.5e5 ;

% ---- boundary cond load ------------
FG = zeros(3*nnod,1) ;
%~ FG( nods2dofs( nodoscargados , 3 ) ) = -1.0 ;
FG = fext* - DIRECTA / (  ( pi * .016 * .5 ) * .07 ) * .2508 * 6.2/15.3 ;
% ===========================================================================================


% ===== 4. Proceso  ==========================================================================

Volumenes = zeros( ntet,1) ;

KG = sparse( 3*nnod, 3*nnod ) ;

for elem = 1:ntet

  nodeselem =  Cellsconec( elem      , 1:4 ) ;
  dofselem  =  nods2dofs ( nodeselem , 3   ) ;
  
  Kml = zeros(12,12) ;

  elecoordmat        = zeros(3,4) ;
  elecoordmat(1,1:4) = nodes( nodeselem , 1 ) ;
  elecoordmat(2,1:4) = nodes( nodeselem , 2 ) ;
  elecoordmat(3,1:4) = nodes( nodeselem , 3 ) ;

  [ deriv , volumen] =  DerivFun( elecoordmat ) ;
  
  if volumen<0, error('Element with negative volume, check connectivity.'), end
  
  Volumenes(elem) = volumen ;
  
  BMat = BMats ( deriv ) ;

  Kml        = BMat' * ConsMat * BMat * volumen ;
  
  KG ( dofselem , dofselem ) = KG ( dofselem , dofselem ) + Kml ;
  
end

Kred = KG(NeumDofs,NeumDofs) ;

Usol = Kred \ FG(NeumDofs)  ;

UG             = zeros( 3*nnod , 1 ) ;
UG( NeumDofs ) = Usol ;



strains  = [];
stresses = [];

for elem = 1:ntet

  nodeselem =  Cellsconec( elem ,1:4) ;

  dofselem = nods2dofs( nodeselem , 3 ) ;

  elecoordmat        = zeros(3,4) ;
  elecoordmat(1,1:4) = nodes( nodeselem , 1 ) ;
  elecoordmat(2,1:4) = nodes( nodeselem , 2 ) ;
  elecoordmat(3,1:4) = nodes( nodeselem , 3 ) ;

  [ deriv , volumen] =  DerivFun( elecoordmat ) ;
  
  BMat = BMats ( deriv ) ;

  strains ( elem , : ) =            ( BMat * UG( dofselem ) )' ;
  stresses( elem , : ) =  ( ConsMat * strains(elem,:)' )' ;
  

end

cellDisp   = cell(2,1);
cellStress = cell(2,1);

cellDisp{1}   = 1 ;
cellDisp{2}   = 'Displacements';
cellDisp{3}   = UG;

cellStress{1} = 2 ;
cellStress{2} = 'Stress';
cellStress{3} = stresses ;

cellDisp
cellStress
mesh2vtk( '../output/ejemploLinearSolid.vtk' , nodes , Cellsconec, cellDisp, cellStress )


function dofs = nods2dofs( nods , degree )
 
  n    = length(nods) ;
  dofs = zeros( n*degree ,1 ) ;

  for i=1:n
    dofs( ((i-1)*degree+1):(i*degree) ) = (nods(i)-1)*degree + (1:degree) ;
  end


% ======================================================================
function [ fun , vol ] = DerivFun( elecoordmat )

  A        = zeros(4,4)   ;
  A(:,1)   = 1.0          ;
  A(:,2:4) = elecoordmat'  ;
  
  invA = inv(A) ;
    
  vol = det(A) / 6.0 ;
  
  fun = invA(2:4,:) ;


% ======================================================================
function BMat = BMats ( deriv )

  BMat = zeros(6,12) ;
  
  for k = 1:4

    for i=1:3
      BMat ( i , (k-1)*3 + i  ) = deriv(i,k) ;
    end

    BMat ( 4 , (k-1)*3 + 2   ) = deriv(3,k)  ;
    BMat ( 4 , (k-1)*3 + 3   ) = deriv(2,k)  ;

    BMat ( 5 , (k-1)*3 + 1   ) = deriv(3,k) ;
    BMat ( 5 , (k-1)*3 + 3   ) = deriv(1,k) ;

    BMat ( 6 , (k-1)*3 + 1   ) = deriv(2,k) ;
    BMat ( 6 , (k-1)*3 + 2   ) = deriv(1,k) ;
      
  end
% ======================================================================




  
function [nodes,Cellsconec, nodoscargados, nodosFijos, fext] = mallayahecha



[ nodes , node_elem , line_elem , trng_elem, tet_elem ] = solid3D_msh_read( '../input/solid3D_maderaUniones.msh' );

fext = zeros(3*size(nodes,1),1) ;

Cellsconec = tet_elem(:,3:end) ;

nodoscargados = [];
nodosFijos    = [];

for i=1:size(trng_elem,1)
  if (trng_elem(i,1)==1001)
    nodosFijos = [ nodosFijos trng_elem(i,2:end) ] ;
  elseif (trng_elem(i,1)==1002)
    nodestrng = trng_elem(i,2:end) ;
    area = 0.5*norm( cross( nodes(nodestrng(2),:) - nodes( nodestrng(1),:) , nodes( nodestrng(3),:) - nodes( nodestrng(1),: ) ) ) ; 

    nodoscargados = [ nodoscargados nodestrng ] ;
    dofs = nods2dofs( nodestrng , 3 ) ;
    fext( dofs(1:3:end) ) = fext( dofs(1:3:end) ) + area * 1/3 ;

  end
end

nodoscargados = unique( nodoscargados ) ;
nodosFijos = unique( nodosFijos);

%~ nodes= [ 0 0 0 4 
%~ 0 0 0.5 4 
%~ 0 0 1 4 
%~ 0.5 0 0 1 
%~ 0.5 0 0.5 1 
%~ 0.5 0 1 1 
%~ 1 0 0 2 
%~ 1 0 0.5 2 
%~ 1 0 1 2 
%~ 0 0.5 0 4 
%~ 0 0.5 0.5 4 
%~ 0 0.5 1 4 
%~ 0.5 0.5 0 0 
%~ 0.5 0.5 0.5 0 
%~ 0.5 0.5 1 0 
%~ 1 0.5 0 2 
%~ 1 0.5 0.5 2 
%~ 1 0.5 1 2 
%~ 0 1 0 4 
%~ 0 1 0.5 4 
%~ 0 1 1 4 
%~ 0.5 1 0 3 
%~ 0.5 1 0.5 3 
%~ 0.5 1 1 3 
%~ 1 1 0 3 
%~ 1 1 0.5 3 
%~ 1 1 1 3 ] ;

%~ nodes(:,4) = [] ;


%~ Cellsconec = [ 2 15 5 14 0 
%~ 2 6 5 15 0 
%~ 2 15 3 6 0 
%~ 1 14 4 13 0 
%~ 1 5 4 14 0 
%~ 1 14 2 5 0 
%~ 2 15 14 11 0 
%~ 2 15 11 12 0 
%~ 2 12 3 15 0 
%~ 1 14 13 10 0 
%~ 1 14 10 11 0 
%~ 1 11 2 14 0 
%~ 5 18 8 17 0 
%~ 5 9 8 18 0 
%~ 5 18 6 9 0 
%~ 4 17 7 16 0 
%~ 4 8 7 17 0 
%~ 4 17 5 8 0 
%~ 5 18 17 14 0 
%~ 5 18 14 15 0 
%~ 5 15 6 18 0 
%~ 4 17 16 13 0 
%~ 4 17 13 14 0 
%~ 4 14 5 17 0 
%~ 11 24 14 23 0 
%~ 11 15 14 24 0 
%~ 11 24 12 15 0 
%~ 10 23 13 22 0 
%~ 10 14 13 23 0 
%~ 10 23 11 14 0 
%~ 11 24 23 20 0 
%~ 11 24 20 21 0 
%~ 11 21 12 24 0 
%~ 10 23 22 19 0 
%~ 10 23 19 20 0 
%~ 10 20 11 23 0 
%~ 14 27 17 26 0 
%~ 14 18 17 27 0 
%~ 14 27 15 18 0 
%~ 13 26 16 25 0 
%~ 13 17 16 26 0 
%~ 13 26 14 17 0 
%~ 14 27 26 23 0 
%~ 14 27 23 24 0 
%~ 14 24 15 27 0 
%~ 13 26 25 22 0 
%~ 13 26 22 23 0 
%~ 13 23 14 26 0  ] ;

%~ Cellsconec(:,5) = [] ;
