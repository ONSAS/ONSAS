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

%function for verification of the numerical results provided by ONSAS, using the analytical expression provided by the user.


function [ numericalVecy ] = analyticSolVerif ...
( analytSol, analyticFunc, loadFactors, controlDisps, timesVec, analyticCheckTolerance, analyticSolFlag, problemName, printflag, outputdir );

  fprintf('----------------------------------------------- \n')
  fprintf('  Analytical solution verification ... ')

  if analyticSolFlag == 1
    analyticalVecy  = analyticFunc(timesVec);
    numericalVecy   = controlDisps ;
    numericalVecx   = timesVec     ;
    analitMagnitude = 'Control displacement' ;
   
  elseif analyticSolFlag == 2
    analyticalVecy  = analyticFunc( controlDisps) ;
    numericalVecy   = loadFactors  ;
    numericalVecx   = controlDisps ;
    analitMagnitude = 'Load Factors' ;
  
  elseif analyticSolFlag == 3
    absError    = ( controlDisps-analytSol ) ;
    normRelativeError = norm( absError ) / norm( analytSol ) ;

  elseif analyticSolFlag == 4
    absError    = ( analytSol - controlDisps ) ;
    normRelativeError = norm( absError  ) / norm( analytSol ) ;

  elseif analyticSolFlag == 5
    absError    = ( controlDisps - analytSol ) ;
    normRelativeError = norm( absError  ) / norm( analytSol ) ;

  end

  if analyticSolFlag == 1 || analyticSolFlag == 2
    nonZeroEntries    = find( analyticalVecy ~= 0 ) ;
    absError          = abs( numericalVecy - analyticalVecy ) ;
    normRelativeError = sum( absError) / sum ( abs(analyticalVecy) ) ; 
  end
  
  % ----------------------------------------
  if analyticSolFlag == 1 || analyticSolFlag == 2 

    nvals = length( numericalVecx); nmaxvistos = 10 ;
    indsMs = (1:nvals)';
    
    if nvals > nmaxvistos
      indsMs = round( linspace(1, nvals, nmaxvistos) ) ;
    end
    
    lw = 2.0 ; ms = 10 ; plotfontsize = 22 ;
    
    figaux = figure ; hold on, grid on
    plot( numericalVecx(indsMs(1)), numericalVecy(indsMs(1))  , 'b-x', 'linewidth', lw,'markersize',ms)
    plot( numericalVecx(indsMs(1)), analyticalVecy(indsMs(1)) , 'r-s', 'linewidth', lw,'markersize',ms)

    plot( numericalVecx(indsMs), numericalVecy( indsMs) , 'bx', 'linewidth', lw,'markersize',ms)
    plot( numericalVecx(indsMs), analyticalVecy(indsMs) , 'rs', 'linewidth', lw,'markersize',ms)

    plot( numericalVecx, numericalVecy  , 'b', 'linewidth', lw,'markersize',ms)
    plot( numericalVecx, analyticalVecy , 'r', 'linewidth', lw,'markersize',ms)
    labx = xlabel('step/time');  laby = ylabel(analitMagnitude) ;
    legend('numeric', 'analytic','location','North')
    set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
    set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;

    currdir = pwd;
    cd(outputdir )
    if printflag == 1
      print( [ problemName '_analyticVerif'  ] ,'-depslatex') ;
    elseif printflag == 2
      print( [ problemName '_analyticVerif.png' ] ,'-dpng') ;
    end
    cd(currdir)
  
    if printflag > 0  
      close(figaux);
    end
  
  end

  if normRelativeError > analyticCheckTolerance ;
    normRelativeError
    error('error: large difference between analytical and numerical solutions!') ;
  else
    fprintf('PASSED.\n')
    fprintf('  Numerical solution error: %12.4e < %10.2e \n',normRelativeError, analyticCheckTolerance)
  end

  fprintf('----------------------------------------------- \n')
end
