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

% function for computation of the normal force and tanget matrix 
% of 3D truss elements using engineering strain.


function [Finte, KTe, stress, dstressdeps, strain ] = ...
  elementTrussInternForce( Xe, Ue, hyperelasparams, A , paramout )

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

  % --- strain ---
  if hyperelasparams(1) == 1
    strain = b1 * Ue ; % small displacements eng. strain

  elseif hyperelasparams(1) == 2
    strain = 0.5 * ( ldef^2 - lini^2 ) / ( lini^2 ) ; % green-lagrange
  
  elseif hyperelasparams(1) == 3
    strain = ( ldef^2 - lini^2 ) / ( lini * (lini + ldef) ) ; % rotated eng
  end
  
  % --- stress and constitutive tensor ---
  [ stress, dstressdeps  ] = hyperElasModels ( strain, hyperelasparams ) ;
  
  TTcl              = Bdif' * e1def ;


  % --- internal forces computation ---
  if hyperelasparams(1) == 1,

    Finte =  A * stress * lini * b1' ;

  elseif hyperelasparams(1) == 2
    
    b2 = 1/(lini^2) * Ue' * Ge ;

    Finte =  A*stress*lini * (b1+b2)' ;
  
  elseif hyperelasparams(1) == 3
    Finte = stress * A * TTcl ;

  end
  % ----------------------------------
    
  % --- tangent matrix computation ---
  if paramout == 2
    
    if hyperelasparams(1) == 1
      KTe   = dstressdeps * A * lini * b1' * b1 ;
 
    elseif hyperelasparams(1) == 2
      KTe   =   stress      * A / lini * Ge  ...
              + dstressdeps * A * lini * ( (b1 + b2)' * (b1 + b2) ) ;

    elseif hyperelasparams(1) == 3
      KMe   = dstressdeps * A / lini * (                TTcl * (TTcl') ) ;
      Ksige =      stress * A / ldef * ( Bdif' * Bdif - TTcl * (TTcl') ) ;
      KTe   = KMe + Ksige ;
    
    end
  else
    KTe = [] ;
  end
  
  Finte = {Finte};
  KTe   = {KTe};
