%~ Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli  

%~ This file is part of ONSAS.

%~ ONSAS is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation, either version 3 of the License, or
%~ (at your option) any later version.

%~ ONSAS is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%~ GNU General Public License for more details.

%~ You should have received a copy of the GNU General Public License
%~ along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.


%Script for generation of the plots of the deformed structure.

function outputPlots( matUts, coordsElemsMat, plotParamsVector, ...
  Conec, Nodes, constantFext, variableFext, strucsize, ...
  controlDisps, visualloadfactor, linearDeformedScaleFactor, ...
  printflag, outputdir, problemName, loadFactors, sectPar, ...
  nonLinearAnalysisBoolean, dynamicAnalysisBoolean, dispsElemsMat, ...
  timeIncr, cellStress, matNts, indexesElems, plotsViewAxis )

	
	if size(matUts,2) == 1 && size(matNts,2) == 1 
		matUts = [zeros(size(matUts,1),1) matUts] ;
		matNts = [zeros(size(matNts,1),1) matNts] ;
	end
		
  nTimesTotal = size( matUts, 2 ) ;

  timeVals = (0:(nTimesTotal-1))*timeIncr;
  
  if (length(plotParamsVector)>1)
    timesPlotsVec = round( linspace(1, nTimesTotal, plotParamsVector(2) ) ) ;
  else
    timesPlotsVec = 1: size(matUts,2) ;
  end

  if plotParamsVector(1) < 3
    outputOctavePlots

  else
    outputVtk

  end  
  
  if length(loadFactors)>1
		if nonLinearAnalysisBoolean || dynamicAnalysisBoolean ~= 0
			outputLoadVSDisplacementsPlot
		end	
  end

end
