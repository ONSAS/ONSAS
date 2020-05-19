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

% ======================================================================

function [systemDeltauRHS, FextG] = computeRHS( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, dispIter, constantFext, variableFext, userLoadsFilename, currLoadFactor, nextLoadFactor, numericalMethodParams, neumdofs, FintGk) 

  [ solutionMethod, stopTolDeltau,   stopTolForces, ...
  stopTolIts,     targetLoadFactr, nLoadSteps,    ...
  incremArcLen, deltaT, deltaNW, AlphaNW, finalTime ] ...
      = extractMethodParams( numericalMethodParams ) ;


  if strcmp( userLoadsFilename , '')
    FextUser = zeros(size(constantFext)) ;
  else
    if solutionMethod == 2
      error('user load not valid for this implementation of arc-length');
    end
    FextUser = feval( userLoadsFilename, nextLoadFactor)  ;
  end


  if solutionMethod == 1

    FextG  = variableFext * nextLoadFactor + constantFext  + FextUser ;

    Resred          = FintGk(neumdofs) - FextG(neumdofs) ;
    systemDeltauRHS = - ( Resred ) ;

  elseif solutionMethod == 2

    FextG  = variableFext * currLoadFactor + constantFext ;

    Resred = FintGk(neumdofs) - FextG(neumdofs)  ;

    % incremental displacement
    systemDeltauRHS = [ -Resred  variableFext(neumdofs) ] ;


  elseif solutionMethod == 3

    [a0NM, a1NM, a2NM, a3NM, a4NM, a5NM, a6NM, a7NM ] = coefsNM( AlphaNW, deltaNW, deltaT ) ;

    FextG  = variableFext * nextLoadFactor + constantFext  + FextUser ;
  
    Fine      =   massMat(neumdofs,neumdofs) * ...
              ( a0NW * ( Uk(neumdofs) - Ut(neumdofs) )  - a2NW * Udott(neumdofs) - a3NW * Udotdott(neumdofs)  )  ;
      
    Fhat      = FextG -Fine ...
                + dampingMat(neumdofs,neumdofs)*...
                (a1NW*(Ut(neumdofs)-Uk(neumdofs)) + Udott(neumdofs)*a4NW  + a5NW*Udotdott(neumdofs))    ...
                - FintGk(neumdofs)                                                                                   ;

  end
    
