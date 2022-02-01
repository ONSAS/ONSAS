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
%
%## Compute RHS
%# this functions computes the right-hand-side of the linear system defined by the numerical method in use.
%#
%#

function [systemDeltauRHS, FextG, fs, Stress, nexTimeLoadFactors ] = computeRHS( modelProperties, BCsData, Ut, Udott, Udotdott, Utp1, Udottp1, Udotdottp1, nextTime, nexTimeLoadFactors )


  [fs, ~, ~ ] = assembler ( ...
    modelProperties.Conec, modelProperties.elements, modelProperties.Nodes, modelProperties.materials, BCsData.KS, Utp1, Udottp1, Udotdottp1, modelProperties.analysisSettings, [1 0 0], modelProperties.nodalDispDamping, nextTime ) ;

  % TO BE REMOVEd!!!
  Stress = [] ;

  Fint = fs{1} ;  Fvis =  fs{2};  Fmas = fs{3} ; Faero = fs{4} ;
  if strcmp( modelProperties.analysisSettings.methodName, 'newtonRaphson' )

    [FextG, nexTimeLoadFactors ]  = computeFext( BCsData.factorLoadsFextCell, BCsData.loadFactorsFuncCell, modelProperties.analysisSettings, nextTime, length(Fint), BCsData.userLoadsFilename, [] ) ;

    systemDeltauRHS = - ( Fint( BCsData.neumDofs ) - FextG( BCsData.neumDofs ) - Faero( BCsData.neumDofs ) ) ;
  elseif strcmp( modelProperties.analysisSettings.methodName, 'arcLength' )

    [FextG, nexTimeLoadFactors ]  = computeFext( BCsData.factorLoadsFextCell, BCsData.loadFactorsFuncCell, modelProperties.analysisSettings, nextTime, length(Fint), BCsData.userLoadsFilename, nexTimeLoadFactors ) ;

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

  elseif strcmp( modelProperties.analysisSettings.methodName, 'newmark' )

    [FextG, nexTimeLoadFactors ]  = computeFext( BCsData.factorLoadsFextCell, BCsData.loadFactorsFuncCell, modelProperties.analysisSettings, nextTime, length(Fint), BCsData.userLoadsFilename, [] ) ;

    rhat      =   Fint ( BCsData.neumDofs ) ...
                + Fvis ( BCsData.neumDofs ) ...
                + Fmas ( BCsData.neumDofs ) ...
                - FextG( BCsData.neumDofs ) ...
                - Faero(BCsData.neumDofs  );

    systemDeltauRHS = -rhat ;

  elseif strcmp( modelProperties.analysisSettings.methodName, 'alphaHHT' )

    fs = assembler ( ...
      modelProperties.Conec, modelProperties.elements, modelProperties.Nodes, modelProperties.materials, BCsData.KS, Ut, Udott, Udotdott, modelProperties.analysisSettings, [1 0 0], modelProperties.nodalDispDamping, nextTime - modelProperties.analysisSettings.deltaT ) ;

    Fintt = fs{1} ;  Fvist =  fs{2};  Fmast = fs{3} ; Faerot = fs{4} ;

    %[FextG, nexTimeLoadFactors ]  = computeFext( BCsData.factorLoadsFextCell, BCsData.loadFactorsFuncCell, modelProperties.analysisSettings, nextTime, length(Fint), BCsData.userLoadsFilename ) ;

    %[FextG, nexTimeLoadFactors ]  = computeFext( BCsData.factorLoadsFextCell, BCsData.loadFactorsFuncCell, modelProperties.analysisSettings, nextTime, length(Fint), BCsData.userLoadsFilename ) ;

%currTime
%nextTime
    [FextG, nexTimeLoadFactors ]  = computeFext( BCsData.factorLoadsFextCell, BCsData.loadFactorsFuncCell, modelProperties.analysisSettings, nextTime, length(Fint), BCsData.userLoadsFilename, [], nextTime ) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
    % fix me
    %%%%%%%%%%%%%%%%%%%%%%%%%%=00000000=============================
    FextGt = FextG ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%=00000000=============================
    %%%%%%%%%%%%%%%%%%%%%%%%%%=00000000=============================
    %%%%%%%%%%%%%%%%%%%%%%%%%%=00000000=============================


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

  end
