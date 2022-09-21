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

function [systemDeltauMatrix, systemDeltauRHS, FextG, fs, nexTimeLoadFactors ] = system_assembler( modelProperties, BCsData, Ut, Udott, Udotdott, Utp1, Udottp1, Udotdottp1, nextTime, nexTimeLoadFactors, previous_state_mat )

  analysisSettings = modelProperties.analysisSettings ;
  nodalDispDamping = modelProperties.nodalDispDamping ;
  neumdofs = BCsData.neumDofs ;
	
  [fs, ~, mats, ~ ] = assembler( modelProperties.Conec, modelProperties.elements, modelProperties.Nodes, modelProperties.materials, BCsData(1).KS, Utp1, Udottp1, Udotdottp1, analysisSettings, [1 0 1 0], nodalDispDamping, nextTime, previous_state_mat  ) ;

  Fint = fs{1} ;  Fvis =  fs{2};  Fmas = fs{3} ; Faero = fs{4} ; 
  
  KT   = mats{1} ; 

  if strcmp( analysisSettings.methodName, 'newmark' ) || strcmp( analysisSettings.methodName, 'alphaHHT' )
    dampingMat = mats{2} ;
    massMat    = mats{3} ;

    global exportFirstMatrices;
    if exportFirstMatrices == true
      KTred      = KT( neumdofs, neumdofs );
      massMatred = massMat(neumdofs,neumdofs);
      save('-mat', 'output/matrices.mat', 'KT','massMat','neumdofs' );
      figure
      spy(full(KT))
      figure
      spy(full(massMat))
      fprintf('matrices exported.\n--------\n')
      exportFirstMatrices = false;
    end
  end


  if strcmp( modelProperties.analysisSettings.methodName, 'newtonRaphson' )

    [FextG, nexTimeLoadFactors ]  = computeFext( modelProperties, BCsData, nextTime, length(Fint), [] ) ;

    systemDeltauRHS = - ( Fint( BCsData.neumDofs ) - FextG( BCsData.neumDofs ) - Faero( BCsData.neumDofs ) ) ;

    systemDeltauMatrix = KT ( neumdofs, neumdofs ) ;
	
  % -----------------------------------------------------------------------------------
  % newmark
  % -----------------------------------------------------------------------------------
  elseif strcmp( modelProperties.analysisSettings.methodName, 'arcLength' )

    [FextG, nexTimeLoadFactors ]  = computeFext( modelProperties, BCsData, nextTime, length(Fint), nexTimeLoadFactors ) ;

    foundLoadCase = false ;
    loadCase = 1 ;
    while ~foundLoadCase
      if isempty( BCsData.factorLoadsFextCell{loadCase} )
        loadCase = loadCase + 1 ;
      else
        foundLoadCase = true ;
      end
    end

    systemDeltauRHS = [ -(Fint(BCsData.neumDofs)-FextG(BCsData.neumDofs)) ...
                        BCsData.factorLoadsFextCell{loadCase}(BCsData.neumDofs) ] ;

    systemDeltauMatrix = KT ( neumdofs, neumdofs ) ;

  % -----------------------------------------------------------------------------------
  % newmark
  % -----------------------------------------------------------------------------------
  elseif strcmp( modelProperties.analysisSettings.methodName, 'newmark' )

    [FextG, nexTimeLoadFactors ]  = computeFext( modelProperties, BCsData, nextTime, length(Fint), [] ) ;

    rhat      =   Fint ( BCsData.neumDofs ) ...
                + Fvis ( BCsData.neumDofs ) ...
                + Fmas ( BCsData.neumDofs ) ...
                - FextG( BCsData.neumDofs ) ...
                - Faero(BCsData.neumDofs  );

    systemDeltauRHS = -rhat ;

    alphaNM = analysisSettings.alphaNM ;
    deltaNM = analysisSettings.deltaNM ;
    deltaT  = analysisSettings.deltaT  ;

    systemDeltauMatrix =                                   KT(         neumdofs, neumdofs ) ...
                         + 1/( alphaNM * deltaT^2)       * massMat(    neumdofs, neumdofs ) ...
                         + deltaNM / ( alphaNM * deltaT) * dampingMat( neumdofs, neumdofs )  ;

  % -----------------------------------------------------------------------------------
  % alpha-HHT
  % -----------------------------------------------------------------------------------
  elseif strcmp( modelProperties.analysisSettings.methodName, 'alphaHHT' )

    if norm( previous_state_mat(:,3) ) >0, error('warning. HHT method with plastic analysis not validated yet'); end
    
    fs = assembler ( ...
      modelProperties.Conec, modelProperties.elements, modelProperties.Nodes, modelProperties.materials, BCsData.KS, Ut, Udott, Udotdott, modelProperties.analysisSettings, [1 0 0 0], modelProperties.nodalDispDamping, nextTime - modelProperties.analysisSettings.deltaT, previous_state_mat  ) ;

    Fintt = fs{1} ;  Fvist =  fs{2};  Fmast = fs{3} ; Faerot = fs{4} ;

    [FextG, nexTimeLoadFactors ]  = computeFext( modelProperties, BCsData, nextTime, length(Fint), [] ) ;

    FextGt = FextG ;

    alphaHHT = modelProperties.analysisSettings.alphaHHT ;

    rhat   =  ( 1 + alphaHHT ) * ( ...
                + Fint  ( BCsData.neumDofs ) ...
                + Fvis  ( BCsData.neumDofs ) ...
                - FextG ( BCsData.neumDofs ) ...
                - Faero ( BCsData.neumDofs ) ...
              ) ...
              ...
              - alphaHHT * ( ...
                + Fintt  ( BCsData.neumDofs ) ...
                + Fvist  ( BCsData.neumDofs ) ...
                - FextGt ( BCsData.neumDofs ) ...
                - Faerot ( BCsData.neumDofs ) ...
                ) ...
              ...
              + Fmas    ( BCsData.neumDofs ) ;

    systemDeltauRHS = -rhat ;


    alphaHHT = analysisSettings.alphaHHT ;
    deltaT   = analysisSettings.deltaT  ;

    deltaNM = (1 - 2 * alphaHHT ) / 2 ;
    alphaNM = (1 - alphaHHT ^ 2 ) / 4 ;

    systemDeltauMatrix = (1 + alphaHHT )                                 * KT         ( neumdofs, neumdofs ) ...
                       + (1 + alphaHHT ) * deltaNM / ( alphaNM*deltaT  ) * dampingMat ( neumdofs, neumdofs )  ...
                       +                         1 / ( alphaNM*deltaT^2) * massMat    ( neumdofs, neumdofs ) ;


  end


