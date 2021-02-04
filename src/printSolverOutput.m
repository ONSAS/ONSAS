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


% Prints analysis output 

function printSolverOutput( outputdir, problemName, timeIndex, lineData )

incrementsResultsFilename = [ outputdir  problemName '_incrementsOutput.tex' ] ;
incrementsNormalForce 		= [ outputdir  problemName '_incrementsNormalForceOutput.tex' ] ;
%~ incrementsTimePerformance = [ outputdir  problemName '_timePerformanceOutput.tex' ] ;

headerIncrements  = [ '$\\#t$ & $ \\lambda(t)$ & its & $\\| RHS \\|$ & $\\| \\Delta u \\|$ & flagExit ' ...
                      ' & npos & nneg  \\\\ \\hline \n \\endhead \n'] ;
%
headerNormalForce = [ 'timeInd & t & $ \\lambda(t)$ & $N_{max}$ & $N_{min}$ \\\\ \\toprule \n'] ;

timeStepEndLine   = [ '\\hdashline\n' ...
                      '%4i & %9.2e & %4i &           &           & %2i & %3i & %3i \\\\ \n' ] ;
%
timeStepIterLine  = [ '     &           & %4i & %9.2e & %9.2e &    &     &     \\\\ \n' ] ;

if timeIndex == 1 && lineData(1)~=1

  % opens and rewrites files
  fileIncrements = fopen( incrementsResultsFilename ,'w');
  fileNormalForce = fopen( incrementsNormalForce, 'w' ) ;

  % write headers
  fprintf( fileIncrements, headerIncrements );
  fprintf(fileNormalForce, headerNormalForce  );
  
else
  fileIncrements 			= fopen( incrementsResultsFilename ,'a' ) ;
  fileNormalForce 		= fopen( incrementsNormalForce, 'a' ) 		;
  %~ fileTimePerformance = fopen( incrementsTimePerformance, 'a' ) ;
end

  %~ fileTimePerformance = fopen( incrementsTimePerformance, 'w' ) ;
  %~ fprintf(fileTimePerformance, [ 'timeInd & t & Solver time (s) & Stores time (s)  \\\\ \\toprule \n'] );
  

%~ timeIndex, currTime, currLoadFactor,  auxIO.itersPerTime, max( max( abs( modelCurrState.Strainst) )*100 ) , ...
  %~ factor_crit , nKeigpos, nKeigneg

%~ stop  
% latex table output
if lineData(1) == 1
  fprintf( fileIncrements, timeStepIterLine, lineData(2), lineData(3), lineData(4) ) ;
  %~ fileIncrements
  %~ lineData
  %~ fclose( fileIncrements )
  %~ timeIndex, currLoadFactor,  auxIO.itersPerTime, max( max( abs( modelCurrState.Strainst) )*100 ) , ...
  %~ factor_crit , nKeigpos, nKeigneg )

%~ fprintf(fileNormalForce, [ ' %4i & %12.3e & %12.3e  & %12.3e & %12.3e \\\\\n' ], ...
  %~ timeIndex, currTime, currLoadFactor, max(currentNormalForces), min(currentNormalForces)  )

%~ fprintf(fileTimePerformance, [' %4i & %12.3e & %5.3e & %5.3e \\\\\n' ], ...
	%~ timeIndex, currTime, tCallSolver, tStores)     
%~ % -----------------------------------

%~ if max( abs( Strainst) ) > 0.05,
  %~ fprintf('WARNING: at timeStep %5i, elements with strain level %4.1f%%!\n', timeIndex, max( abs( Strainst) )*100 ),
  
elseif lineData(1) == 2

%~ printSolverOutput( outputDir, problemName, timeIndex, [ 2 nextLoadFactor dispIter deltaErrLoad norm(deltaured) stopCritPar nKeigpos nKeigneg ] ) ;
  fprintf( fileIncrements, timeStepEndLine, timeIndex,  lineData(2), lineData(3), lineData(4), lineData(5), lineData(6) );
  
end

% close files
fclose(fileIncrements);
fclose(fileNormalForce);
%~ fclose(fileTimePerformance);

%~ if ( stopTimeIncrBoolean == 1 )
  %~ fprintf('----------------------------------------------- \n');
%~ end
