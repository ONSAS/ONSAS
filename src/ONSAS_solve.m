% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, J. Bruno Bazzano,
% Joaquin Viera, Marcelo Forets, Jean-Marc Battini. 
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
function [ matUs, loadFactorsMat, matFinte ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData )
%md
%md init structures to store solutions
matUs          = modelCurrSol.U              ;
loadFactorsMat = modelCurrSol.currLoadFactorsVals ;
matUdots       = modelCurrSol.Udot           ;
cellStress     = { modelCurrSol.Stress }     ;
matFinte 			 = modelCurrSol.Finte 				 ;
%md
%md#### Incremental time analysis
%md sets stopping boolean to false
finalTimeReachedBoolean = false ;
%mdand starts the iteration
fprintf('Starting time analysis. Time index: ')
while finalTimeReachedBoolean == false

 
  %if mod(modelCurrSol.timeIndex,10)==0,
    fprintf(' %3i,', modelCurrSol.timeIndex),
  %end

  % compute the model state at next time
  modelNextSol = timeStepIteration( modelCurrSol, modelProperties, BCsData ) ;

  % check if final time was reached
  finalTimeReachedBoolean = ( modelNextSol.currTime - modelProperties.analysisSettings.finalTime ) ...
                        >= ( -(modelProperties.analysisSettings.finalTime) * 1e-8 ) ;

  % store results and update structs
  modelCurrSol   =   modelNextSol ;
  matUs          = [ matUs          modelCurrSol.U                   ] ;
  loadFactorsMat = [ loadFactorsMat ; modelCurrSol.currLoadFactorsVals ] ;
	matFinte 			 = [ matFinte modelCurrSol.Finte ] ; 	
	
  % generate vtk file for the new state
  if strcmp( modelProperties.plotsFormat, 'vtk' )
    vtkMainWriter( modelCurrSol, modelProperties );
  end % if vtk output format

end %while time

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
