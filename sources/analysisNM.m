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


% This function implements the Newmark's method for the analysis of one time step. 
% The input includes the structural properties, the displacements at time t
% and the loads for the next time tp1. 
% The output includes the displacements, velocities, accelerations and internal forces.

function ...
%  outputs ---
[ Utp1, Udottp1, Udotdottp1, FintGtp1, dispIter, Strainsk, Stressk ] ...
  = analysisNM ( ...
  % --- inputs ---
  % constant data
  Conec, secGeomProps, coordsElemsMat, neumdofs, nnodes, hyperElasParamsMat, ...
  constantFext, variableFext, KS, ...  %
  massMat, dampingMat, a0NW, ...
  a1NW, a2NW, a3NW, a4NW,...
   a5NW, a6NW,   a7NW, ...
  % model variable data
  Ut, Udott, Udotdott, nextLoadFactor, stopTolDeltau, stopTolForces, stopTolIts, userLoadsFilename, nextTime ) ;
% -------------------------------------------------------


% nelems and dofs per node
nelems    = size(Conec,1) ;
ndofpnode = 6             ;

% ---  initialize parameters ---
iterDispConverged = 0                         ;
dispIter          = 0                         ;
deltau            = zeros(length(neumdofs),1) ;

if strcmp( userLoadsFilename ,'')
  FextUser = zeros(size(constantFext)) ;
else
  FextUser = feval( userLoadsFilename, nextTime)  ;
end

Uk          = Ut ;
[FintGt]    = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Ut,[], 1  ) ;%  el 1 es por la fuerza
FextGtp1    = variableFext * nextLoadFactor + constantFext  + FextUser                                     ;
FextGtp1red = FextGtp1(neumdofs)                                                                           ;                                                          

Fine      =  massMat(neumdofs,neumdofs) * ...
              ( a0NW * ( Uk(neumdofs) - Ut(neumdofs) )  - a2NW * Udott(neumdofs) - a3NW * Udotdott(neumdofs)  ) ;
              
Fhat        = FextGtp1red  -  Fine ...
              + dampingMat(neumdofs,neumdofs)*...
              (a1NW*(Ut(neumdofs)-Uk(neumdofs)) + Udott(neumdofs)*a4NW  + a5NW*Udotdott(neumdofs))    ...
              - FintGt(neumdofs)          ;                                                                       
       

while ( iterDispConverged == 0 )
  dispIter += 1;
  
  % computes tangent matrix
  [~, KT]  = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk,[], 2  ) ;
  
  %VALSKT= abs(diag (eig(KT)))

  %computes deltau and rfresh Ut
  Kroof       = KT(neumdofs, neumdofs) + (massMat*a0NW + a1NW*dampingMat) (neumdofs, neumdofs) ;                     
  
  %VALSKROOF= abs(diag (eig(Kroof)) )
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  deltaured   = Kroof \ (Fhat);                                                                                     
  %%%%%%%%%%%%%%%%%%%%%%%%

  normadeltau = norm(deltaured);                                                                                    
  Uk(neumdofs)= Uk(neumdofs ) + deltaured ;
  normaUk     = norm( Uk(neumdofs ) ) ;                                                                             

  % updates model variables and computes internal, hat and inertial forces
  [FintGk, ~, Strainsk, Stressk ]  = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk,[], 1  ); %1 para sacar la fuerza
  
  Fine      =   massMat(neumdofs,neumdofs) * ...
              ( a0NW * ( Uk(neumdofs) - Ut(neumdofs) )  - a2NW * Udott(neumdofs) - a3NW * Udotdott(neumdofs)  )  ;

              
  Fhat      = FextGtp1red -Fine ...
              + dampingMat(neumdofs,neumdofs)*...
              (a1NW*(Ut(neumdofs)-Uk(neumdofs)) + Udott(neumdofs)*a4NW  + a5NW*Udotdott(neumdofs))    ...
              - FintGk(neumdofs)                                                                                   ;

  % --- stopping criteria verification ---
  logicDispStop = ( normadeltau < ( normaUk * stopTolDeltau ) ) ;
  logicForcStop = ( norm(Fhat)  < ( norm(FextGtp1(neumdofs) )  * stopTolForces ) ) ;
                
  if logicForcStop
    stopCritPar = 1 ;      iterDispConverged = 1 ;

  elseif logicDispStop
    stopCritPar = 2 ;      iterDispConverged = 1 ;

  elseif ( dispIter >= stopTolIts )
    warning('displacements iteration stopped by max iterations.');
    stopCritPar = 3 ;      iterDispConverged = 1 ;
  end
  
  fprintf(' %12.3e   %03i  %12g\n', norm(Fhat), dispIter, nextTime )
  % -------------------------

end

Utp1       = Uk                                         ;
Udotdottp1 = a0NW*(Utp1-Ut) - a2NW*Udott - a3NW*Udotdott;
Udottp1    = Udott + a6NW*Udotdott + a7NW*Udotdottp1    ;
FintGtp1   = FintGk                                     ;
