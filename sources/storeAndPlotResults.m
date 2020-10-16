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

deltaT    = modelNextSol.currTime - modelCurrSol.currTime ;
timeIndex = modelCurrSol.timeIndex ; 

if timeIndex == 1
  matUs      = modelCurrSol.U          ;
  matUdots   = modelCurrSol.Udot       ;
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

if exist( 'iniMatUs' )
  if timeIndex+2 <= size( iniMatUs,2)
    Utp10 = iniMatUs(:,timeIndex+2) ;
  else
    Utp10 = [] ;
  end
end

indsNormal = [ find(elementsParamsMat(Conec(:,4+2)) == 2 ) ; find( elementsParamsMat(Conec(:,4+2)) == 3 ) ]' ;


if timeIndex == 1
    normalForcesIni = zeros( nElems, 1 ) ;
    sigxs = modelCurrSol.Stress(:,1) ;
    if ~isempty(indsNormal) == 1
        for indexArea = 1:(size( crossSecsParamsMat,1))
            typeSec          = crossSecsParamsMat(indexArea,1)                  ;
            indexElemTypeSec = find(elementsParamsMat(Conec(:,4+4)) == typeSec) ;
            typeSecParams    = crossSecsParamsMat(indexArea,:)                  ;
            if typeSecParams(1) == 1 %general section
                areaTypeSec  = typeSecParams( 2 )                                ;
            elseif typeSecParams(1) == 2 %rectangular section
                areaTypeSec  = typeSecParams(2)*typeSecParams(3)                 ;
            elseif typeSecParams(1) == 3 %circular section
                diameter     = typeSecParams(2)                                  ;
                areaTypeSec  = pi*diameter^2/4                                   ;
            else
                error(' section type not implemented yet, please create an issue')
            end
            normalForcesIni( indexElemTypeSec ) =  sigxs(indexElemTypeSec) .* areaTypeSec ;
        end
    end
end
% initiate normalForces with initial values
normalForces = normalForcesIni ;

if ~isempty(indsNormal) == 1
    for indexArea = 1:(size(crossSecsParamsMat,1))                    
        typeSec          = crossSecsParamsMat(indexArea,1)                    ;
        indexElemTypeSec = find(elementsParamsMat(Conec(:,4+4)) == typeSec)     ;
        typeSecParams    = crossSecsParamsMat(indexArea,:)                      ;
        if typeSecParams(1) == 1 %general section
            areaTypeSec = typeSecParams( 2 ) ;
        elseif typeSecParams(1) == 2 %rectangular section
            areaTypeSec = typeSecParams(2)*typeSecParams(3)      ;
        elseif typeSecParams(1) == 3 %circular section
            diameter = typeSecParams(2) ;
            areaTypeSec = pi*diameter^2/4           ;
        else
            error(' section type not implemented yet, please create an issue')
        end
        normalForces( indexElemTypeSec ) =  sigxs(indexElemTypeSec) .* areaTypeSec ;
    end
end



if length(indsNormal) > 0
    sigxs = modelNextSol.Stress(:,1) ;
    areasVector       = zeros(nElems,1);
        for indexElem = 1:nElems ;
            typeSec = Conec(indexElem,4+4)                        ;
            elemCrossSecParams = crossSecsParamsMat(typeSec,:)   ;
            if elemCrossSecParams(1) == 1 %general section
                areasVector(indexElem) = elemCrossSecParams( 2 ) ;
            elseif elemCrossSecParams(1) == 2 %rectangular section
                areasVector(indexElem) = elemCrossSecParams(2)*elemCrossSecParams(3)      ;
            elseif elemCrossSecParams(1) == 3
                diameter = elemCrossSecParams(2) ;
                areasVector(indexElem) = pi*diameter^2/4           ;
            else
                error(' section type not implemented yet, please create an issue')
            end
        end
  normalForces( indsNormal ) =  sigxs .* areasVector ;
end

if (storeBoolean == 1)
  matUs      (:, timeIndex+1 )       = modelNextSol.U                  ;
  matUdots   (:, timeIndex+1 )       = modelNextSol.Udot               ;
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
    [ vtkNodes, vtkDispMat, vtkNormalForces, vtkStress ...
      , vtkConec, elem2VTKCellMap ] ...
      = vtkGeometry( ...
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
