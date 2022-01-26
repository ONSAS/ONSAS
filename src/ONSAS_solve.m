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
fprintf('Starting time analysis. Time index: ')
while finalTimeReachedBoolean == false

  %if mod(modelCurrSol.timeIndex,10)==0,
    fprintf(' %3i,', modelCurrSol.timeIndex),
  %end

  % compute the model state at next time
  modelNextSol = timeStepIteration( modelCurrSol, modelProperties, BCsData ) ;

  % check if final time was reached
  finalTimeReachedBoolean = ( modelNextSol.currTime - modelProperties.analysisSettings.finalTime ) ...
                        >= ( -(modelProperties.analysisSettings.finalTime) * 1e-8 ) ;

  % store results and update structs
  modelCurrSol   =   modelNextSol ;
  matUs          = [ matUs          modelCurrSol.U                   ] ;
  loadFactorsMat = [ loadFactorsMat ; modelCurrSol.currLoadFactorsVals ] ;

modelCurrSol
stop
  % generate vtk file for the new state
  if strcmp( modelProperties.plotsFormat, 'vtk' )
    vtkMainWriter( modelCurrSol, modelProperties );
  end % if vtk output format

end %while time
%md
