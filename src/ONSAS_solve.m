% Copyright (C) 2020, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
%   Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini, Sebastian Toro
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

%md### ONSAS_solve
%md Function that performs the time analysis with the model structs as input.
%md
function [ matUs, loadFactorsMat ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData )
%md
%md init structures to store solutions
matUs          = modelCurrSol.U              ;
loadFactorsMat = modelCurrSol.currLoadFactorsVals ;
matUdots       = modelCurrSol.Udot           ;
cellStress     = { modelCurrSol.Stress }     ;
%md
%md#### Incremental time analysis
%md sets stopping boolean to false
finalTimeReachedBoolean = false ;
%mdand starts the iteration
while finalTimeReachedBoolean == false

  % compute the model state at next time
  modelNextSol = timeStepIteration( modelCurrSol, modelProperties, BCsData ) ;

  % check if final time was reached
  finalTimeReachedBoolean = ( modelNextSol.currTime - modelProperties.analysisSettings.finalTime ) ...
                        >= ( -(modelProperties.analysisSettings.finalTime) * 1e-8 ) ;

  % store results and update structs
  modelCurrSol   =   modelNextSol ;
  matUs          = [ matUs          modelCurrSol.U                   ] ;
  loadFactorsMat = [ loadFactorsMat ; modelCurrSol.currLoadFactorsVals ] ;

  % generate vtk file for the new state
  if ( modelProperties.plotParamsVector(1) == 3  ) && ...
     ( sum( find( modelCurrSol.timeIndex == timesPlotsVec ) ) == 1

    vtkMainWriter( modelCurrSol, )
  end

end %while time
%md
