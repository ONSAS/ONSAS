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


% Prints analysis output 

incrementsResultsFilename = [ outputdir  problemName '_incrementsOutput.tex' ] ;
incrementsNormalForce 		= [ outputdir  problemName '_incrementsNormalForceOutput.tex' ] ;
incrementsTimePerformance 	= [ outputdir  problemName '_timePerformanceOutput.tex' ] ;


if currTime == 0
  fprintf('----------------------------------------------- \n')
  % prints header table line
  fileIncrements = fopen( incrementsResultsFilename ,'w') ; 
  fprintf(fileIncrements, [ 'timeInd & t & $ \\lambda(t)$ & dispits & maxStrain (\\%%) & BucklingFac ' ...
  ' & npos & nneg  \\\\ \\toprule \n'] );
  fileNormalForce = fopen( incrementsNormalForce, 'w' ) ;
  fprintf(fileNormalForce, [ 'timeInd & t & $ \\lambda(t)$ & $N_{max}$ & $N_{min}$ \\\\ \\toprule \n'] );
  fileTimePerformance = fopen( incrementsTimePerformance, 'w' ) ;
  fprintf(fileTimePerformance, [ 'timeInd & t & Solver time (s) & Stores time (s)  \\\\ \\toprule \n'] );
else
  fileIncrements 			= fopen( incrementsResultsFilename ,'a' ) ;
  fileNormalForce 		= fopen( incrementsNormalForce, 'a' ) 		;
  fileTimePerformance = fopen( incrementsTimePerformance, 'a' ) ;
end

% latex table output
fprintf(fileIncrements, [ ' %4i & %12.3e & %12.3e  & %4i  & %5.2f & %12.5e & %5i & %3i \\\\\n' ], ...
  timeIndex, currTime, currLoadFactor,  auxIO.itersPerTime, max( abs( modelCurrState.Strainst) )*100 , ...
  factor_crit , nKeigpos, nKeigneg )
fprintf(fileNormalForce, [ ' %4i & %12.3e & %12.3e  & %12.3e & %12.3e \\\\\n' ], ...
  timeIndex, currTime, currLoadFactor, max(currentNormalForces), min(currentNormalForces)  )
fprintf(fileTimePerformance, [' %4i & %12.3e & %5.3e & %5.3e \\\\\n' ], ...
	timeIndex, currTime, tCallSolver, tStores)     
% -----------------------------------

if max( abs( Strainst) ) > 0.05,
  fprintf('WARNING: at timeStep %5i, elements with strain level %4.1f%%!\n', timeIndex, max( abs( Strainst) )*100 ),
end

fclose(fileIncrements);
fclose(fileNormalForce);
fclose(fileTimePerformance);



if ( stopTimeIncrBoolean == 1 )
  %~ totaliterationstimeinseconds = toc
  fprintf('----------------------------------------------- \n');
end

