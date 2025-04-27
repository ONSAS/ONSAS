% Copyright 2025, ONSAS Authors (see documentation)
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
function [booleanConverged, stopCritPar, deltaErrLoad, normFext] = convergenceTest( ...
                                                                                   analysisSettings, redFext, redDeltaU, redUk, dispIter, systemDeltauRHS)

  stopTolDeltau = analysisSettings.stopTolDeltau;
  stopTolForces = analysisSettings.stopTolForces;
  stopTolIts    = analysisSettings.stopTolIts;

  normaUk       = norm(redUk);
  normadeltau   = norm(redDeltaU);

  deltaErrLoad  = norm(systemDeltauRHS);
  normFext      = norm(redFext);

  logicDispStop = (normadeltau  < (normaUk  * stopTolDeltau));
  logicForcStop = (deltaErrLoad < ((normFext + (normFext < stopTolForces)) * stopTolForces))  * (deltaErrLoad > 0);

  if isnan(norm(redDeltaU))
    error(' ERROR: the displacements are NOT A NUMBER!');
  end

  if logicForcStop
    stopCritPar = 1;
    booleanConverged = 1;
  elseif logicDispStop
    stopCritPar = 2;
    booleanConverged = 1;
  elseif dispIter >= stopTolIts
    stopCritPar = 3;
    booleanConverged = 1;
  else
    stopCritPar = 0;
    booleanConverged = 0;
  end
