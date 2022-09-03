% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% Joaquin Viera, J. Bruno Bazzano, Marcelo Forets, Jean-Marc Battini.
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

%md### ONSAS_solve
%md Function that performs the time analysis with the model structs as input.
%md
function [ matUs, loadFactorsMat, cellFint ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData )
%md
%md init structures to store solutions
matUs          = modelCurrSol.U              ;
loadFactorsMat = modelCurrSol.currLoadFactorsVals ;
matUdots       = modelCurrSol.Udot           ;
cellStress     = { modelCurrSol.Stress }     ;

cellFint = {};


%md
%md#### Incremental time analysis
%md sets stopping boolean to false
finalTimeReachedBoolean = false ;
%mdand starts the iteration
fprintf('|                                                 |\n')
fprintf('| Analysis progress:   |0       50       100| %%   |\n')
fprintf('|                      |')
plotted_bars = 0 ;

plots_counter = 0 ;

iterations_average = 0 ;
iterations_maximum = 0 ;
iterations_strop_crit_vec = [ 0 0 0 ] ;

tic
while finalTimeReachedBoolean == false

  percent_time = round( (modelCurrSol.timeIndex*modelProperties.analysisSettings.deltaT) ...
                       / modelProperties.analysisSettings.finalTime * 20 ) ;

  while plotted_bars < percent_time,
    fprintf('=')
    plotted_bars = plotted_bars +1 ;
  end

  % compute the model state at next time
  modelNextSol = timeStepIteration( modelCurrSol, modelProperties, BCsData ) ;

  % iterations statistics
  iterations_average = ...
    (   iterations_average* (modelNextSol.timeIndex-2) + modelNextSol.timeStepIters ) ...
    / ( modelNextSol.timeIndex-1) ;

  iterations_maximum = max( iterations_maximum, modelNextSol.timeStepIters ) ;
  iterations_strop_crit_vec( modelNextSol.timeStepStopCrit ) = ...
    iterations_strop_crit_vec( modelNextSol.timeStepStopCrit ) + 1 ;
  % ------------------------------



  % check if final time was reached
  finalTimeReachedBoolean = ( modelNextSol.currTime - modelProperties.analysisSettings.finalTime ) ...
                        >= ( -(modelProperties.analysisSettings.finalTime) * 1e-8 ) ;


 
  % store results and update structs
  modelCurrSol   	=  	modelNextSol ;
  matUs          	= [ matUs          modelCurrSol.U                   ] ;
  loadFactorsMat 	= [ loadFactorsMat ; modelCurrSol.currLoadFactorsVals ] ;
	
	if length(cellFint) == 0
		cellFint{1} = zeros(size(modelCurrSol.matFint)) ;
	end
		
	cellFint{end+1}	= modelCurrSol.matFint ;
	
		
  % generate vtk file for the new state
  if strcmp( modelProperties.plots_format, 'vtk' )
    vtkMainWriter( modelCurrSol, modelProperties );
  end % if vtk output format

end %while time
time_solve = toc ;
fprintf('|     |\n')


% ---- print iteration statistics -----
fprintf('| Time: %6.1f sec                                |\n',time_solve)
fprintf('|                                                 |\n')
fprintf('| Iters:  avg  max | Stop by: force  disp  iters  |\n')
fprintf('|        %4.1f  %3i |          %5i %5i  %5i  |\n', ...
  iterations_average, iterations_maximum, iterations_strop_crit_vec(1), ...
  iterations_strop_crit_vec(2), iterations_strop_crit_vec(3) )
% -------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BLOQUE DE ANALISIS MODAL PROVISORIO %%%%%%
global modalAnalysisBoolean
if ~isempty(modalAnalysisBoolean) && modalAnalysisBoolean
  pwd
  genpath( [ pwd '/output'])
  addpath( genpath( [ pwd '/output'] ) ); load( 'matrices.mat' ) ;
  Kred = KT(BCsData.neumDofs,BCsData.neumDofs);
  Mred = massMat(BCsData.neumDofs,BCsData.neumDofs);
  Mred = Mred + speye(size(Mred,1));
  numModes = 10;
  [PHI, OMEGA] = eigs(Mred^(-1)*Kred,numModes,'sm');

  modelPropertiesModal = modelProperties ;
  modelCurrSolModal    = modelCurrSol    ;

  for i = 1:4
    fprintf(' generating mode %2i vtk\n', i) ;
    modelPropertiesModal.problemName = [ modelProperties.problemName sprintf('_mode_%02i_', i ) ] ;
    modelCurrSolModal.U = zeros( size(modelCurrSol.U, 1) , 1 )    ;
    modelCurrSolModal.U( BCsData.neumDofs ) = PHI(:,i)  ;
    vtkMainWriter( modelCurrSolModal, modelPropertiesModal ) ;
  end

  save('-binary','Modal.mat','PHI','OMEGA')
  fprintf(' MODAL ANALYSIS DONE. Setting modalAnalysisBoolean to false.\n')
  modalAnalysisBoolean = false ;

end %endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%md
