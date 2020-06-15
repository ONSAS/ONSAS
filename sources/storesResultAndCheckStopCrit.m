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


% script for updating and storing variables at each time increment.
% In this script the dispsElemesMat is also created, with the displacements 
% of all nodes, including rotation of releases. 
% In this script the analysis stopping criteria is checked.

% --------------   computes magnitudes for next step  -------------------------- 
Utp1 = modelNextState.Ut ;
Ut   = modelCurrState.Ut ;

Udott    = modelNextState.Udott    ;
Udotdott = modelNextState.Udotdott ;
 
% updates displacements
convDeltau = Utp1 - Ut ;

timeIndex = modelCurrState.timeIndex ;

loadFactors  ( timeIndex +1 )  = BCsData.nextLoadFactor  ;
controlDisps ( timeIndex +1 )  = Utp1( controlDofsAndFactors(1) ) * controlDofsAndFactors(2) ;
timesVec     ( timeIndex +1 )  = deltaT * timeIndex ;
% ------------------------------------------------------------------------------


% ----------------   updates data structures and time --------------------------
if dynamicAnalysisBoolean == 0
%~ BCsData.nextLoadFactor
%~ stop
  currLoadFactor = BCsData.nextLoadFactor    ;
  nextLoadFactor = currLoadFactor + numericalMethodParams(5) / nLoadSteps ;
else

  currLoadFactor = loadFactorsFunc( modelCurrState.currTime + deltaT ) ;
  nextLoadFactor = loadFactorsFunc( modelCurrState.currTime + 2*deltaT ) ;
  
end

BCsData.currLoadFactor = currLoadFactor ;
BCsData.nextLoadFactor = nextLoadFactor ;

modelCurrState            = modelNextState ; 
modelCurrState.convDeltau = convDeltau     ;

%~ modelCurrState.Udott 
%~ modelCurrState.Udotdott 

%~ timeIndex      = timeIndex +1      ;
%~ currTime       = currTime + deltaT ;

% ------------------------------------------------------------------------------

  %~ targetLoadFactr = loadFactorsFunc( finalTime ) ;
  %~ nLoadSteps      = round( finalTime / deltaT ) ;
  %~ factor_crit     = 0;
  %~ nKeigpos        = 0;
  %~ nKeigneg        = 0 ;



% ---------------------------------------------------


% --- stores displacements ---
matUts = [ matUts modelNextState.Ut ] ;

tangentMatricesCell{timeIndex} = modelNextState.systemDeltauMatrix ;
 
% normal forces calculation

%~ indselems12 = find( ( Conec(:,7) == 1) | ( Conec(:,7) == 2) ) ;
Areas = secGeomProps(Conec(:,6),1) ;
currentNormalForces = modelCurrState.Stresst(:,1) .* Areas ;

matNts = [ matNts currentNormalForces ] ;

% ---------------------------------------------------
% prints step analysis data/results

itersPerTimeVec( timeIndex )    = auxIO.itersPerTime ;

if dynamicAnalysisBoolean == 0
  factor_crit = modelNextState.factorCrit ;
  nKeigneg = modelNextState.nKeigneg ;
  nKeigpos= modelNextState.nKeigpos;
else
  factor_crit = 0 ;
  nKeigneg = 0 ;
  nKeigpos = 0 ;
end

% ---------------       evals stop time incr crit          ---------------------
%~ if dynamicAnalysisBoolean == 1
%~ modelNextState.currTime 
%~ finalTime

  if ( modelNextState.currTime > finalTime )
    stopTimeIncrBoolean = 1 ; fprintf('%4i.\n',timeIndex);
  else
    while contProgr < ( timeIndex / ( nLoadSteps*.05 ) )
      contProgr = contProgr + 1 ;
      fprintf('=')
    end
  end
%~ else
  %~ if ( (nextLoadFactor - numericalMethodParams(5)) > 1e-6*numericalMethodParams(5) ) || ( timeIndex > nLoadSteps ) % || ( abs( currTime - finalTime) < (deltaT*1e-4) )
    %~ stopTimeIncrBoolean = 1 ; fprintf('%4i.\n',timeIndex);
  %~ end
%~ end
% ------------------------------------------------------------------------------

%~ stop

