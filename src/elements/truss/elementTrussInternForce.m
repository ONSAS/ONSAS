% Copyright 2023, ONSAS Authors (see documentation)
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
% 
% function for computation of the normal force and tanget matrix
% of 3D truss elements using engineering strain.


function [Finte, KTe, stress, dstressdeps, strain, acum_plas_strain ] = ...
  elementTrussInternForce( Xe, Ue, modelName, modelParams, A, previous_state )
  
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

  acum_plas_strain =  previous_state{3} ;

  % --- strain ---
  if strcmp( modelName, 'elastic-linear')
    strain = b1 * Ue ; % small displacements

  elseif strcmp( modelName, 'elastic-SVK')
    strain = 0.5 * ( ldef^2 - lini^2 ) / ( lini^2 ) ;  % green-lagrange

  elseif strcmp( modelName, 'elastic-rotEngStr') || strcmp( modelName, 'plastic-rotEngStr')
    strain = ( ldef^2 - lini^2 ) / ( lini * (lini + ldef) ) ; % rotated eng

  elseif strcmp( modelName, 'elastic-rotLogStr') || strcmp( modelName, 'plastic-rotLogStr')
    strain = log(ldef/lini) ; % logarithmic strain

  else
    modelName
    error('strain expression not known for given material modelName');
  end

  % --- stress ---
  if strcmp( modelName, 'elastic-linear')

    E           = modelParams(1) ;
    stress      = E * strain ;
    Finte       = A * stress * lini * b1' ;
    KTe         = E * A * lini * b1' * b1 ;
    dstressdeps = E ;
  
  elseif strcmp( modelName, 'elastic-SVK')
  
    lambda = modelParams(1) ;
    mu     = modelParams(2) ;
  
    E           = mu*(3*lambda+2*mu)/(lambda+mu) ;
    stress      = E * strain ;
    dstressdeps = E ;
  
    b2    = 1/(lini^2) * Ue' * Ge ;
    Finte =  A*stress*lini * (b1+b2)' ;
  
    KTe   =   stress      * A / lini * Ge  ...
            + dstressdeps * A * lini * ( (b1 + b2)' * (b1 + b2) ) ;
  
  elseif strcmp( modelName, 'elastic-rotEngStr') || strcmp( modelName, 'elastic-rotLogStr')
  
    E           = modelParams(1) ;
  
    stress      = E * strain ;
    dstressdeps = E ;
       
    Finte = stress * A * TTcl ;
  
    KMe   = dstressdeps * A / lini * (                TTcl * (TTcl') ) ;
    Ksige =      stress * A / ldef * ( Bdif' * Bdif - TTcl * (TTcl') ) ;
    KTe   = KMe + Ksige ;  
  
  elseif strcmp( modelName, 'plastic-rotEngStr') || strcmp( modelName, 'plastic-rotLogStr')
      
    E           = modelParams(1) ;
  
    stress_n           = previous_state{1,:}(1)  ;
    strain_n           = previous_state{2,:}(1)  ;
    acum_plas_strain_n =  previous_state{3} ;
        
    Kplas       = modelParams(2) ;
    sigma_Y_0   = modelParams(3) ;
    
    stress_Elas = stress_n + E * (strain - strain_n) ;
    phi_tr = abs( stress_Elas ) - ( sigma_Y_0 + Kplas*acum_plas_strain_n ) ;
    
    if phi_tr < 0 % elastic behavior
      stress           = stress_Elas ;
      dstressdeps      = E ;
      acum_plas_strain = acum_plas_strain_n ; 
          
    else % elasto-plastic behavior
      delta_gamma = phi_tr / ( E + Kplas ) ;
  
      stress           = stress_Elas - E*delta_gamma * sign( stress_Elas ) ;
      dstressdeps      = E*Kplas / ( E + Kplas ) ;
      acum_plas_strain = acum_plas_strain_n + delta_gamma ; 
    end
  
    Finte = stress * A * TTcl ;
  
    KMe   = dstressdeps * A / lini * (                TTcl * (TTcl') ) ;
    Ksige =      stress * A / ldef * ( Bdif' * Bdif - TTcl * (TTcl') ) ;
    KTe   = KMe + Ksige ;
    
  end

  % ----------------------------------

  Finte = {Finte};
  KTe   = {KTe};
