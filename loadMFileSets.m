% ======================================================================
% Actualize GUI values
% ======================================================================



% Analysis options
set(handles.numericalTitle, 'visible', 'on'), set(handles.numMethod, 'visible', 'on')

if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0

	buttonString = 'Linear Analysis' ;
	set(handles.numMethod, 'value', 1)
	set(handles.linAna, 'value', 1)
	if exist('numericalMethodParams')~=0
		handles.numVec = 1 ;
	end

elseif dynamicAnalysisBoolean == 1
	buttonString = 'Dynamic NonLinear Analysis' ;
	set(handles.dynNonLinAna, 'value', 1)
elseif nonLinearAnalysisBoolean == 1 
	buttonString = 'NonLinear Analysis' ;
	set(handles.nonLinAna, 'value', 1) ;

end

if strcmp(buttonString, 'Linear Analysis')
	set(handles.nonLinAna, 'value', 0), set(handles.dynNonLinAna, 'value', 0)
	handles.anaVec = [ 1 ] ;
	set(handles.scaleText, 'visible', 'on'), set(handles.scale, 'visible', 'on')
	set(handles.controlDofText, 'visible', 'off'), set(handles.controlDof, 'visible', 'off')
    
	set(handles.numMethod, 'value', 1)
	set(handles.numMethod, 'string', {'Choose method', 'User Function'} )
	set(handles.tolDeltaUText, 'visible', 'off'), set(handles.tolDeltaU, 'visible', 'off') 
	set(handles.tolForcesText, 'visible', 'off'), set(handles.tolForces, 'visible', 'off') 
	set(handles.tolIterText, 'visible', 'off'), set(handles.tolIter, 'visible', 'off')
	set(handles.targetLoadFacText, 'visible', 'off'), set(handles.targetLoadFac, 'visible', 'off') 
	set(handles.loadStepsText, 'visible', 'off'), set(handles.loadSteps, 'visible', 'off')
	set(handles.incremALText, 'visible', 'off'), set(handles.incremAL, 'visible', 'off')
	
	set(handles.timeIncrText, 'visible', 'off'), set(handles.timeIncr, 'visible', 'off')
	set(handles.finalTimeText, 'visible', 'off'), set(handles.finalTime, 'visible', 'off')
	set(handles.tolDeltaUDynText, 'visible', 'off'), set(handles.tolDeltaUDyn, 'visible', 'off')
	set(handles.tolForcesDynText, 'visible', 'off'), set(handles.tolForcesDyn, 'visible', 'off')
	set(handles.tolIterDynText, 'visible', 'off'), set(handles.tolIterDyn, 'visible', 'off')
	set(handles.deltaNWText, 'visible', 'off'), set(handles.deltaNW, 'visible', 'off')
	set(handles.alphaNWText, 'visible', 'off'), set(handles.alphaNW, 'visible', 'off')	
    
elseif strcmp(buttonString, 'Dynamic NonLinear Analysis')
	set(handles.linAna, 'value', 0), set(handles.nonLinAna, 'value', 0)
	
	handles.anaVec = [ 2 ] ;
	set(handles.scaleText, 'visible', 'off'), set(handles.scale, 'visible', 'off')
	set(handles.controlDofText, 'visible', 'on'), set(handles.controlDof, 'visible', 'on')
	
	set(handles.numMethod, 'string', {'Choose method', 'Newmark'} )
	set(handles.numMethod, 'value', 2)
	set(handles.deltaTText, 'visible', 'off'), set(handles.deltaT, 'visible', 'off')
	set(handles.finalTimeLinText, 'visible', 'off'), set(handles.finalTimeLin, 'visible', 'off')
	set(handles.indexTimeSolText, 'visible', 'off'), set(handles.indexTimeSol, 'visible', 'off')
	
	set(handles.tolDeltaUText, 'visible', 'off'), set(handles.tolDeltaU, 'visible', 'off') 
	set(handles.tolForcesText, 'visible', 'off'), set(handles.tolForces, 'visible', 'off') 
	set(handles.tolIterText, 'visible', 'off'), set(handles.tolIter, 'visible', 'off')
	set(handles.targetLoadFacText, 'visible', 'off'), set(handles.targetLoadFac, 'visible', 'off') 
	set(handles.loadStepsText, 'visible', 'off'), set(handles.loadSteps, 'visible', 'off')
	set(handles.incremALText, 'visible', 'off'), set(handles.incremAL, 'visible', 'off')
	if get(handles.numMethod, 'value') == 2
		set(handles.timeIncrText, 'visible', 'on'), set(handles.timeIncr, 'visible', 'on')
		set(handles.finalTimeText, 'visible', 'on'), set(handles.finalTime, 'visible', 'on')
		set(handles.tolDeltaUDynText, 'visible', 'on'), set(handles.tolDeltaUDyn, 'visible', 'on')
		set(handles.tolForcesDynText, 'visible', 'on'), set(handles.tolForcesDyn, 'visible', 'on')
		set(handles.tolIterDynText, 'visible', 'on'), set(handles.tolIterDyn, 'visible', 'on')
		set(handles.deltaNWText, 'visible', 'on'), set(handles.deltaNW, 'visible', 'on')
		set(handles.alphaNWText, 'visible', 'on'), set(handles.alphaNW, 'visible', 'on')
		handles.numVec = 1 ;
	else
		set(handles.timeIncrText, 'visible', 'off'), set(handles.timeIncr, 'visible', 'off')
		set(handles.finalTimeText, 'visible', 'off'), set(handles.finalTime, 'visible', 'off')
		set(handles.tolDeltaUDynText, 'visible', 'off'), set(handles.tolDeltaUDyn, 'visible', 'off')
		set(handles.tolForcesDynText, 'visible', 'off'), set(handles.tolForcesDyn, 'visible', 'off')
		set(handles.tolIterDynText, 'visible', 'off'), set(handles.tolIterDyn, 'visible', 'off')
		set(handles.deltaNWText, 'visible', 'off'), set(handles.deltaNW, 'visible', 'off')
		set(handles.alphaNWText, 'visible', 'off'), set(handles.alphaNW, 'visible', 'off')
		handles.numVec = 0 ;
	end
elseif strcmp(buttonString, 'NonLinear Analysis')  
	set(handles.linAna, 'value', 0) , set(handles.dynNonLinAna, 'value', 0) 
	handles.anaVec = [ 3 ] ;
	set(handles.scaleText, 'visible', 'off'), set(handles.scale, 'visible', 'off')
	set(handles.controlDofText, 'visible', 'on'), set(handles.controlDof, 'visible', 'on')
	
	
	set(handles.numMethod, 'string', {'Choose method', 'Newton Raphson', 'Newton Raphson-Arc Length'} )
	
	set(handles.deltaTText, 'visible', 'off'), set(handles.deltaT, 'visible', 'off')
	set(handles.finalTimeLinText, 'visible', 'off'), set(handles.finalTimeLin, 'visible', 'off')
	set(handles.indexTimeSolText, 'visible', 'off'), set(handles.indexTimeSol, 'visible', 'off')
	set(handles.timeIncrText, 'visible', 'off'), set(handles.timeIncr, 'visible', 'off')
	set(handles.finalTimeText, 'visible', 'off'), set(handles.finalTime, 'visible', 'off')
	set(handles.tolDeltaUDynText, 'visible', 'off'), set(handles.tolDeltaUDyn, 'visible', 'off')
	set(handles.tolForcesDynText, 'visible', 'off'), set(handles.tolForcesDyn, 'visible', 'off')
	set(handles.tolIterDynText, 'visible', 'off'), set(handles.tolIterDyn, 'visible', 'off')
	set(handles.deltaNWText, 'visible', 'off'), set(handles.deltaNW, 'visible', 'off')
	set(handles.alphaNWText, 'visible', 'off'), set(handles.alphaNW, 'visible', 'off')  
	
	if get(handles.numMethod, 'value') == 2
		set(handles.tolDeltaUText, 'visible', 'on'), set(handles.tolDeltaU, 'visible', 'on') 
		set(handles.tolForcesText, 'visible', 'on'), set(handles.tolForces, 'visible', 'on') 
		set(handles.tolIterText, 'visible', 'on'), set(handles.tolIter, 'visible', 'on')
		set(handles.targetLoadFacText, 'visible', 'on'), set(handles.targetLoadFac, 'visible', 'on') 
		set(handles.loadStepsText, 'visible', 'on'), set(handles.loadSteps, 'visible', 'on')
		set(handles.incremALText, 'visible', 'off'), set(handles.incremAL, 'visible', 'off') 
		handles.numVec = 1 ;
	elseif get(handles.numMethod, 'value') == 3
		set(handles.incremALText, 'visible', 'on'), set(handles.incremAL, 'visible', 'on')
		handles.numVec = 2 ;
	else
		set(handles.tolDeltaUText, 'visible', 'off'), set(handles.tolDeltaU, 'visible', 'off') 
		set(handles.tolForcesText, 'visible', 'off'), set(handles.tolForces, 'visible', 'off') 
		set(handles.tolIterText, 'visible', 'off'), set(handles.tolIter, 'visible', 'off')
		set(handles.targetLoadFacText, 'visible', 'off'), set(handles.targetLoadFac, 'visible', 'off') 
		set(handles.loadStepsText, 'visible', 'off'), set(handles.loadSteps, 'visible', 'off')
		set(handles.incremALText, 'visible', 'off'), set(handles.incremAL, 'visible', 'off')
		handles.numVec = 0 ;
	end
	
end


%~ handles.prevValTolDeltaU = stopTolDeltau ;
%~ set(handles.tolDeltaU, 'string', num2str(handles.prevValTolDeltaU))
