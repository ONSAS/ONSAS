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


%script for updating and storing variables at each time increment. In this script the dispsElemesMat is also created, with the displacements of all nodes, including rotation of releases. In this script the analysis stopping criteria is checked.

if dynamicAnalysisBoolean == 0
  deltaT    = targetLoadFactr/nLoadSteps ;
  finalTime = targetLoadFactr ;
end

% --------------   computes magnitudes for next step  -------------------------- 
Utp1 = modelNextState.Ut ;
Ut   = modelCurrState.Ut ;
 
% updates displacements
convDeltau = Utp1 - Ut ;

loadFactors  ( timeIndex +1 )  = BCsNextState.nextLoadFactor  ;
controlDisps ( timeIndex +1 )  = Utp1(controlDof)*controlDofFactor ;
timesVec     ( timeIndex +1 )  = deltaT * timeIndex ;
% ------------------------------------------------------------------------------


% ----------------   updates data structures and time --------------------------
if dynamicAnalysisBoolean == 0
  currLoadFactor = BCsNextState.nextLoadFactor    ;
  nextLoadFactor = currLoadFactor + targetLoadFactr / nLoadSteps ;
end

BCsNextState.currLoadFactor = currLoadFactor ;
BCsNextState.nextLoadFactor = nextLoadFactor ;

modelCurrState            = modelNextState ; 
modelCurrState.convDeltau = convDeltau     ;

if timeIndex == 1,
  fprintf('Time/step: %4i, ',timeIndex);
elseif mod( timeIndex, 20) == 0,
  fprintf('%4i,',timeIndex);
end

timeIndex      = timeIndex +1      ;
currTime       = currTime + deltaT ;

% ------------------------------------------------------------------------------

  %~ targetLoadFactr = loadFactorsFunc( finalTime ) ;
  %~ nLoadSteps      = round( finalTime / deltaT ) ;
  %~ factor_crit     = 0;
  %~ nKeigpos        = 0;
  %~ nKeigneg        = 0 ;



% ---------------------------------------------------


% updates load factor

% stores displacements
matUts = [ matUts modelNextState.Ut ] ;

% normal forces calculation

indselems12 = find( ( Conec(:,7) == 1) || ( Conec(:,7) == 2) ) ;
Areas = secGeomProps(Conec(:,6),1) ;
currentNormalForces = modelCurrState.Stresst(:) .* Areas ;

matNts = [ matNts currentNormalForces ] ;

% ---------------------------------------------------
% prints step analysis data/results

itersPerTimeVec( timeIndex )    = auxIO.itersPerTime ;

if dynamicAnalysisBoolean == 0
  factor_crit = modelCurrState.factorCrit ;
  nKeigneg = modelCurrState.nKeigneg ;
  nKeigpos= modelCurrState.nKeigpos;
else
  factor_crit = 0 ;
  nKeigneg = 0 ;
  nKeigpos = 0 ;
end

printsOutputScreen


% ---------------       evals stop time incr crit          ---------------------
if dynamicAnalysisBoolean == 1
  if ( timeIndex > nLoadSteps ) || ( abs( currTime - finalTime) < (deltaT*1e-10))
    stopTimeIncrBoolean = 1 ; fprintf('%4i.\n',timeIndex);
  end
else
  if (nextLoadFactor > targetLoadFactr) || ( timeIndex > nLoadSteps ) % || ( abs( currTime - finalTime) < (deltaT*1e-4) )
    stopTimeIncrBoolean = 1 ; fprintf('%4i.\n',timeIndex);
  end
end
% ------------------------------------------------------------------------------
