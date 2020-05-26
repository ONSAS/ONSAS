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


function [ booleanConverged, stopCritPar, deltaErrLoad ] = convergenceTest( ...
  numericalMethodParams, redFint, redFext, redDeltaU, redUk, dispIter ) 

  [ solutionMethod, stopTolDeltau,   stopTolForces, ...
    stopTolIts,     targetLoadFactr, nLoadSteps,    ...
    incremArcLen, deltaT, deltaNW, AlphaNW, finalTime ] ...
        = extractMethodParams( numericalMethodParams ) ;

  normaUk       = norm( redUk )               ;
  normadeltau   = norm( redDeltaU         )   ;
  deltaErrLoad  = norm( redFint - redFext )   ;
  normFext      = norm( redFext )             ;
  
  logicDispStop = ( normadeltau  < ( normaUk  * stopTolDeltau ) )  ;
  logicForcStop = ( deltaErrLoad < ( normFext * stopTolForces ) )  ;
                
  if logicForcStop
    stopCritPar = 1 ;      booleanConverged = 1 ;
    warning('stopped by forcesssss')

  elseif logicDispStop
    stopCritPar = 2 ;      booleanConverged = 1 ;
    fprintf('\nstopped by displacements\n\n')

  elseif ( dispIter >= stopTolIts )
    warning('displacements iteration stopped by max iterations.');
    stopCritPar = 3 ;      booleanConverged = 1 ;
  else
    booleanConverged = 0;  stopCritPar = [];
  end
  
