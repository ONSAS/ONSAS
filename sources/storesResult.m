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


deltaT    = modelNextSol.currTime - modelCurrSol.currTime ;
timeIndex = modelCurrSol.timeIndex ; 

if timeIndex == 1
  matUs      = modelCurrSol.U ;
  cellStress = { modelCurrSol.Stress } ;
  
  controlDisps               = 0 ;
  controlDisps(timeIndex, :) = modelCurrSol.U( controlDofsAndFactors(:,1) ) ...
                                            .* controlDofsAndFactors(:,2) ;

  loadFactors                = BCsData.currLoadFactor  ;
end

% ----------------   updates data structures and time --------------------------

loadFactors  ( timeIndex+1 )       = BCsData.nextLoadFactor  ;
controlDisps ( timeIndex+1 )       = modelNextSol.U ( controlDofsAndFactors(:,1) ) ...
                                                   .* controlDofsAndFactors(:,2) ;
timesVec     ( timeIndex+1 )       = deltaT * timeIndex ;


if storeBoolean = 1
  matUs      (:, timeIndex+1 )       = modelNextSol.U                  ;
  cellStress   { timeIndex+1 }       = modelNextSol.Stress             ;
  tangentMatricesCell{ timeIndex+1 } = modelNextSol.systemDeltauMatrix ;
end
% ------------------------------------------------------------------------------

% stores iterations
itersPerTimeVec( timeIndex )    = modelNextSol.timeStepIters ;

while contProgr < ( modelCurrSol.currTime / ( finalTime*.05 ) )
  contProgr = contProgr + 1 ;
  fprintf('=')
end
