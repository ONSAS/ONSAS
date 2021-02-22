
function ONSAS_solve( modelCurrSol, modelProperties, BCsData)

% --- increment step analysis ---
stopTimeIncrBoolean = 0 ;
while ( stopTimeIncrBoolean == 0 )

  % -----   computes the model state at the next load/time step   -----
  [ modelNextSol, BCsData ] = timeStepIteration ( modelCurrSol, BCsData, modelProperties );

  % --- checks stopping criteria and stores model state

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
            typeSec          = crossSecsParamsMat(indexArea,1)   ;
            indexElemTypeSec = find((Conec(:,4+4)) == indexArea) ;
            typeSecParams    = crossSecsParamsMat(indexArea,:)   ;
            if typeSecParams(1) == 1 %general section
                areaTypeSec  = typeSecParams( 2 )                ;
            elseif typeSecParams(1) == 2 %rectangular section
                areaTypeSec  = typeSecParams(2)*typeSecParams(3) ;
            elseif typeSecParams(1) == 3 %circular section
                diameter     = typeSecParams(2)                  ;
                areaTypeSec  = pi*diameter^2/4                   ;
            else
                error(' section type not implemented yet, please create an issue')
            end
            normalForcesIni( indexElemTypeSec ) =  sigxs(indexElemTypeSec) .* areaTypeSec ;
        end
    end
end
% initiate normalForces with initial values
normalForces = normalForcesIni ;
sigxs = modelCurrSol.Stress(:,1) ;
if ~isempty(indsNormal) == 1
    for indexArea = 1:(size(crossSecsParamsMat,1))                    
        typeSec          = crossSecsParamsMat(indexArea,1)   ;
        indexElemTypeSec = find((Conec(:,4+4)) == indexArea) ;
        typeSecParams    = crossSecsParamsMat(indexArea,:)   ;
        if typeSecParams(1) == 1 %general section
            areaTypeSec = typeSecParams( 2 ) ;
        elseif typeSecParams(1) == 2 %rectangular section
            areaTypeSec = typeSecParams(2)*typeSecParams(3)   ;
        elseif typeSecParams(1) == 3 %circular section
            diameter = typeSecParams(2)                       ;   
            areaTypeSec = pi*diameter^2/4                     ;   
        else
            error(' section type not implemented yet, please create an issue')
        end
        normalForces( indexElemTypeSec ) =  sigxs(indexElemTypeSec) .* areaTypeSec ;
    end
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



  % --- update state data --- 
  modelCurrSol           = modelNextSol ;
  
  BCsData.currLoadFactor = BCsData.nextLoadFactor                   ;
  BCsData.nextLoadFactor = loadFactorsFunc( modelCurrSol.currTime + deltaT ) ;

  % ----   evals stop time incr crit    --------
  stopTimeIncrBoolean = modelCurrSol.currTime >= finalTime ;
    
end

fprintf( '\n| end time-step%4i - (max,avg) iters: (%3i,%5.2f) | \n ',...
  modelCurrSol.timeIndex, max( itersPerTimeVec ) , mean( itersPerTimeVec(2:end) ) );




