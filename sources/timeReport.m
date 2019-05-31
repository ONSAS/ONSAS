% Scripts that prints the time report

fprintf('\n')
fprintf('============================================\n')
fprintf('            SYSTEM TIME REPORT \n')
fprintf('============================================\n\n')

% Reading and vars verif/defs
fprintf('         Reading and vars verif/defs \n')
fprintf('============================================\n')
fprintf('-------------------------------------------------------------\n')
fprintf(' Reading input file:                             |  %5.3fs  |\n', tReadingInput )
fprintf('-------------------------------------------------------------\n')
fprintf(' Variables verification:                         |  %5.3fs  |\n', tVarVer)
fprintf('-------------------------------------------------------------\n')
fprintf(' Input aux defs:                                 |  %5.3fs  |\n', tInputAuxDefs)
fprintf('-------------------------------------------------------------\n')
fprintf('-------------------------------------------------------------\n')
fprintf(' Total elapsed time in reading and verification: |  %5.3fs  |\n', tReadingInput+tVarVer+tInputAuxDefs)
fprintf('-------------------------------------------------------------\n\n')

% Analysis time
fprintf('               Analysis time \n')
fprintf('============================================\n')
if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
	% Geometry reading
	fprintf('-------------------------------------------------------------\n')
	fprintf(' Geometry reading time:                          |  %5.3fs  |\n', tGeomReading)
	fprintf('-------------------------------------------------------------\n')
	% Stiffness matrix assembly
	fprintf(' Stiffness matrix assembly:                      |  %5.3fs  |\n', tStiffMatrix)
	fprintf('-------------------------------------------------------------\n')
	% Loads assembly
	fprintf(' Loads assembly:                                 |  %5.3fs  |\n', tLoadsAssembly)
	fprintf('-------------------------------------------------------------\n')
	% System resolution
	fprintf(' System resolution:                              |  %5.3fs  |\n', tSystemResolution)
	fprintf('-------------------------------------------------------------\n')
	% Elems solic and disps
	fprintf(' Elems solic and disps:                          |  %5.3fs  |\n', tSolicDisps)
	fprintf('-------------------------------------------------------------\n')
	% Total time
	fprintf('-------------------------------------------------------------\n')
	fprintf(' Total elapsed time in linear analysis:          |  %5.3fs  |\n', tGeomReading+tStiffMatrix+tLoadsAssembly+tSystemResolution+tSolicDisps)
	fprintf('-------------------------------------------------------------\n\n')
else
	
end

% Plot time
if plotParamsVector(1)>0

	fprintf('               Plots time \n')
	fprintf('============================================\n')
	fprintf(' Auxiliar defs:                                  |  %5.3fs  |\n', tPrevDefs)
	fprintf('-------------------------------------------------------------\n')
	if plotParamsVector(1) < 3
		fprintf(' Margin defs:                                    |  %5.3fs  |\n', tMarginDef)
		fprintf('-------------------------------------------------------------\n')
		fprintf(' Deformed shape:                                 |  %5.3fs  |\n', tDefShape)
		fprintf('-------------------------------------------------------------\n')
		if size(matUts,2) > 1
			fprintf(' Load factor vs control disp:                    |  %5.3fs  |\n', tLoadFac)
			fprintf('-------------------------------------------------------------\n')
		end
		fprintf(' Normal force:                                   |  %5.3fs  |\n', tNormalForce)
		fprintf('-------------------------------------------------------------\n')
	else
		fprintf('VTK ... todo')
	end
	if length(loadFactors)>1	
		fprintf(' Load vs Dips:                                   |  %5.3fs  |\n', tLoadDisps)
		fprintf('-------------------------------------------------------------\n')
	end
end

% Report time
if reportBoolean == 1
	fprintf('               Report time \n')
	fprintf('============================================\n')
	fprintf(' Report writing:                                 |  %5.3fs  |\n', tReport)
	fprintf('-------------------------------------------------------------\n')
end
