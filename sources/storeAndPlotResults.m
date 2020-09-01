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

  if length( controlDofsAndFactors ) > 0  
    controlDisps               = 0 ;
    controlDisps(timeIndex, :) = modelCurrSol.U( controlDofsAndFactors(:,1) ) ...
                                              .* controlDofsAndFactors(:,2) ;
  end
  
  loadFactors                = BCsData.currLoadFactor  ;
end

% ----------------   updates data structures and time --------------------------

loadFactors  ( timeIndex+1 )       = BCsData.nextLoadFactor  ;
if length( controlDofsAndFactors ) > 0  
  controlDisps ( timeIndex+1 )       = modelNextSol.U ( controlDofsAndFactors(:,1) ) ...
                                                   .* controlDofsAndFactors(:,2) ;
end
timesVec     ( timeIndex+1 )       = deltaT * timeIndex ;

normalForces = zeros( nElems, 1 ) ;
indsNormal = [ find( Conec(:,4+2) == 2 ) ; find( Conec(:,4+2) == 3 ) ]' ;

if timeIndex == 1
  normalForcesIni = zeros( nElems, 1 ) ;
  sigxs = modelCurrSol.Stress(:,1) ;
  if length(indsNormal) > 0
    crossSecParamsMat = cell2mat(crossSecsParams( Conec( indsNormal, 4+4) , 1 ));
    areasVector       = zeros(size(crossSecParamsMat,1),1);
    for indexElem = 1: size(areasVector,1);
        if crossSecParamsMat (indexElem,1)==3
            areasVector(indexElem) = pi*crossSecParamsMat(indexElem,2)^2/4;
        elseif crossSecParamsMat (indexElem,1)==2
             areasVector(indexElem)= crossSecParamsMat(indexElem,2)*crossSecParamsMat(indexElem,3);
        elseif crossSecsParamsElem crossSecParamsMat(indexElem,1)==1
            areasVector(indexElem) = crossSecParamsMat (indexElem,2);
        end
    end
    normalForcesIni( indsNormal ) =  sigxs .* areasVector ;
  end
end

if length(indsNormal) > 0
  sigxs = modelNextSol.Stress(:,1) ;
      crossSecParamsMat = cell2mat(crossSecsParams( Conec( indsNormal, 4+4) , 1 ));
    areasVector       = zeros(size(crossSecParamsMat,1),1);
    for indexElem = 1: size(areasVector,1);
        if crossSecParamsMat (indexElem,1)==3
            areasVector(indexElem) = pi*crossSecParamsMat(indexElem,2)^2/4;
        elseif crossSecParamsMat (indexElem,1)==2
             areasVector(indexElem)= crossSecParamsMat(indexElem,2)*crossSecParamsMat(indexElem,3);
        elseif crossSecsParamsElem crossSecParamsMat(indexElem,1)==1
            areasVector(indexElem) = crossSecParamsMat (indexElem,2);
        end
    end
  normalForces( indsNormal ) =  sigxs .* areasVector ;
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
    [ vtkNodes, vtkDispMat, vtkNormalForces, vtkStress, vtkConec, elem2VTKCellMap ] = vtkGeometry( ...
      modelProperties.coordsElemsMat , Conec, crossSecsParamsMat, modelCurrSol.U, normalForcesIni, Nodes, elementsParamsMat, modelCurrSol.Stress ) ;

    % data
    [ cellPointData, cellCellData, filename ] = vtkData( outputDir, problemName, 1, vtkNormalForces, vtkStress, vtkDispMat ) ;

    % writes file
    vtkWriter( filename, vtkNodes, vtkConec , cellPointData, cellCellData ) ;
  end

  aux = timesPlotsVec == modelNextSol.timeIndex ;
  
  if sum( aux ) > 0

    indplot = find( aux ) ;

    % generates deformed nodes coords
    [ vtkNodes, vtkDispMat, vtkNormalForces, vtkStress ] = vtkGeometry( ...
      modelProperties.coordsElemsMat , Conec, crossSecsParamsMat, modelNextSol.U, normalForces, Nodes, elementsParamsMat, modelNextSol.Stress  ) ;

    % data
    [ cellPointData, cellCellData, filename ] = vtkData( outputDir, problemName, indplot, vtkNormalForces, vtkStress, vtkDispMat ) ;
    
    % writes file
    vtkWriter( filename, vtkNodes, vtkConec , cellPointData, cellCellData ) ;
      
  end
end
