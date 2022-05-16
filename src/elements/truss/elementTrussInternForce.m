% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, J. Bruno Bazzano,
% Joaquin Viera, Marcelo Forets, Jean-Marc Battini. 
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
 
% function for computation of the normal force and tanget matrix
% of 3D truss elements using engineering strain.


function [Finte, KTe, stress, dstressdeps, strain ] = ...
  elementTrussInternForce( Xe, Ue, hyperElasModel, hyperElasParams, A )

  Xe    = Xe'     ;
  Xedef = Xe + Ue ;

  Bdif = [ -eye(3) eye(3) ] ;
  Ge   = Bdif' * Bdif       ;

  % initial/deformed lengths
  lini = sqrt( sum( ( Bdif * Xe    ).^2 ) ) ;
  ldef = sqrt( sum( ( Bdif * Xedef ).^2 ) ) ;

  % normalized reference and deformed co-rotational vector
  e1ref = Bdif * Xe    / lini ;
  e1def = Bdif * Xedef / ldef ;

  b1 = 1/(lini^2) * Xe' * Ge ;

  TTcl              = Bdif' * e1def ;

  % --- strain ---
  if strcmp( hyperElasModel, 'linearElastic')

    strain = b1 * Ue ; % small displacements eng. strain

    E      = hyperElasParams(1) ;
    stress = E * strain ;
    Finte  =  A * stress * lini * b1' ;
    KTe    = E * A * lini * b1' * b1 ;

  elseif strcmp( hyperElasModel, 'SVK')
    strain = 0.5 * ( ldef^2 - lini^2 ) / ( lini^2 ) ;  % green-lagrange

    lambda = hyperElasParams(1) ;
    mu     = hyperElasParams(2)     ;

    E           = mu*(3*lambda+2*mu)/(lambda+mu) ;
    stress      = E * strain ;
    dstressdeps = E ;

    b2    = 1/(lini^2) * Ue' * Ge ;
    Finte =  A*stress*lini * (b1+b2)' ;

    KTe   =   stress      * A / lini * Ge  ...
            + dstressdeps * A * lini * ( (b1 + b2)' * (b1 + b2) ) ;

  elseif strcmp( hyperElasModel, '1DrotEngStrain')
    strain = ( ldef^2 - lini^2 ) / ( lini * (lini + ldef) ) ; % rotated eng

    E           = hyperElasParams(1) ;
    stress      = E * strain ;
    dstressdeps = E ;

    Finte = stress * A * TTcl ;

    KMe   = dstressdeps * A / lini * (                TTcl * (TTcl') ) ;
    Ksige =      stress * A / ldef * ( Bdif' * Bdif - TTcl * (TTcl') ) ;
    KTe   = KMe + Ksige ;

  end


  % ----------------------------------

  Finte = {Finte};
  KTe   = {KTe};
