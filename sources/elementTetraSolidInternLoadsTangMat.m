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

% function for computation of nodal forces and tangent stiffness matrix for 3D 4 nodes tetraedron element

function [ Finte, KTe, strain, stress ] = elementTetraSolidInternLoadsTangMat( tetcoordmat, Ue, params )
  
  Finte = zeros(24,1) ;
  KTe = sparse(24,24) ;
  
  [ deriv , vol ] =  DerivFun( tetcoordmat ) ;
  
  if vol<0, error('Element with negative volume, check connectivity.'), end
  tetVol = vol ;
  BMat = BMats ( deriv ) ;
  
  E  = params(1) ;
  nu = params(2) ;

  mu    = E / (2.0 * ( 1.0 + nu ) ) ;
  Bulk  = E / (3.0 * ( 1.0 - 2.0 * nu ) ) ;
  
  eyetres      = eye(3)     ;
  eyevoig      = zeros(6,1) ;
  eyevoig(1:3) = 1.0        ;

  ConsMat = zeros(6,6) ;
  ConsMat (1:3,1:3) = Bulk  + 2* mu * ( eyetres - 1.0/3.0 ) ;
  ConsMat (4:6,4:6) =            mu *   eyetres ;
    
  Kml    = BMat' * ConsMat * BMat * tetVol ;
  
  KTe([1:2:end], [1:2:end])  =  Kml ;

  strain = BMat * Ue ;
  stress = ConsMat * strain ;
  Fint   = stress' * BMat * tetVol ;
  
  Finte(1:2:end) = Fint ;
  


