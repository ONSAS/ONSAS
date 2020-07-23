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

normalForces = zeros( nElems, 1 ) ;
indsNormal = [ find( Conec(:,7) == 1 ) ; find( Conec(:,7) == 2 ) ]' ;

if timeIndex == 1
  normalForcesIni = zeros( nElems, 1 ) ;
  sigxs = modelCurrSol.Stress(:,1) ;
  if length(indsNormal) > 0
    normalForcesIni( indsNormal ) =  sigxs .* crossSecsParams( Conec( indsNormal, 6) , 1 ) ;
  end
end

if length(indsNormal) > 0
  sigxs = modelNextSol.Stress(:,1) ;
  normalForces( indsNormal ) =  sigxs .* crossSecsParams( Conec( indsNormal, 6) , 1 ) ;
end

if (storeBoolean == 1)
  matUs      (:, timeIndex+1 )       = modelNextSol.U                  ;
  cellStress {   timeIndex+1 }       = modelNextSol.Stress             ;
  tangentMatricesCell{ timeIndex+1 } = modelNextSol.systemDeltauMatrix ;
  matNs      (:, timeIndex+1 )       = normalForces                    ;
end

% stores iterations
itersPerTimeVec( timeIndex+1 )  = modelNextSol.timeStepIters ;

% ------------------------------------------------------------------------------


%~ contProgr
%~ modelNextSol.currTime
%~ ( finalTime*.05 )

while contProgr <= ( modelNextSol.currTime / ( finalTime*.05 ) )
  contProgr = contProgr + 1 ;
  fprintf('=')
end


if plotParamsVector(1) == 3

  if modelCurrSol.timeIndex == 1,
    
    % generate connectivity
    [ vtkNodes, vtkDispMat, vtkNormalForces, vtkConec, elem2VTKCellMap ] = vtkGeometry( ...
      modelProperties.coordsElemsMat , Conec, sectPar, modelCurrSol.U, normalForcesIni, Nodes ) ;

    % data
    [ cellPointData, cellCellData, filename ] = vtkData( outputDir, problemName, 1, vtkNormalForces, {}, vtkDispMat ) ;

    % writes file
    vtkWriter( filename, vtkNodes, vtkConec , cellPointData, cellCellData ) ;
  end

  aux = timesPlotsVec == modelNextSol.timeIndex ;
  
  if sum( aux ) > 0

    indplot = find( aux ) ;

    % generates deformed nodes coords
    [ vtkNodes, vtkDispMat, vtkNormalForces ] = vtkGeometry( ...
      modelProperties.coordsElemsMat , Conec, sectPar, modelNextSol.U, normalForces, Nodes ) ;

    % data
    [ cellPointData, cellCellData, filename ] = vtkData( outputDir, problemName, indplot, vtkNormalForces, {}, vtkDispMat ) ;
    
    % writes file
    vtkWriter( filename, vtkNodes, vtkConec , cellPointData, cellCellData ) ;
      
  end
end
