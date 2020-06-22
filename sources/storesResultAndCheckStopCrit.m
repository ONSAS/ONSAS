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


deltaT = modelNextSol.currTime - modelCurrSol.currTime ;

% --- update state data --- 
modelCurrSol   = modelNextSol ;
BCsData.currLoadFactor = BCsData.nextLoadFactor                            ;
BCsData.nextLoadFactor = loadFactorsFunc( modelCurrSol.currTime + deltaT ) ;


% ----------------   updates data structures and time --------------------------
timeIndex                   = modelCurrSol.timeIndex ;

loadFactors  ( timeIndex )       = BCsData.currLoadFactor  ;
controlDisps ( timeIndex )       = modelCurrSol.U ( controlDofsAndFactors(:,1) ) ...
                                                 .* controlDofsAndFactors(:,2) ;
timesVec     ( timeIndex )       = deltaT * timeIndex ;
matUs                            = [ matUs modelCurrSol.U ] ;
tangentMatricesCell{ timeIndex } = modelCurrSol.systemDeltauMatrix ;

cellStress{ timeIndex } = modelCurrSol.Stress ;

% ------------------------------------------------------------------------------

% stores iterations
itersPerTimeVec( timeIndex )    = modelCurrSol.timeStepIters ;

while contProgr < ( timeIndex / ( nLoadSteps*.05 ) )
  contProgr = contProgr + 1 ;
  fprintf('=')
end
