% Copyright (C) 2021, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
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

function printSolverOutput( outputdir, problemName, lineData )

iterationResultsFile = [ outputdir  problemName '_iterations.tex' ] ;

headerIncrements  = [ '$\\#t$ & $t$ & its & $\\| RHS \\|$ & $\\| \\Delta u \\|$ & flagExit \\\\ \\hline \n \\endhead \n'] ;
%
timeStepIterLine  = [ '     &           & %4i & %9.2e & %9.2e &      \\\\ \n' ] ;
%
timeStepEndLine   = [ '%4i & %9.2e & %4i &           &           & %3i  \\\\ \n \\hdashline \n' ] ;

if lineData(1)==0 % print header

  % opens and rewrites files
  fileIncrements = fopen( iterationResultsFile ,'w');

  % write headers
  fprintf( fileIncrements, headerIncrements ) ;

else
  fileIncrements 			= fopen( iterationResultsFile ,'a' ) ;
end

if lineData(1) == 1 % iteration information
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


elseif lineData(1) == 2 %end of iteration information
  fprintf( fileIncrements, timeStepEndLine, lineData(2), lineData(3), lineData(4), lineData(5) );

end

% close file
fclose(fileIncrements);
