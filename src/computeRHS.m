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

% ======================================================================

function [systemDeltauRHS, FextG, fs, Stress] = computeRHS( modelProperties, BCsData, Ut, Udott, Udotdott, Utp1, Udottp1, Udotdottp1, nextTime )

  fs = assembler ( ...
    modelProperties.Conec, modelProperties.elements, modelProperties.Nodes, modelProperties.materials, BCsData.KS, Utp1, 1, Udottp1, Udotdottp1, modelProperties.analysisSettings ) ;
  
  Fint = fs{1} ;  Fvis =  fs{2};  Fmas = fs{3} ;  

  if strcmp( modelProperties.analysisSettings.methodName, 'newtonRaphson' )

    FextG           = computeFext( BCsData, modelProperties.analysisSettings, nextTime, length(Fint) ) ;

    systemDeltauRHS = - ( Fint( BCsData.neumDofs ) - FextG( BCsData.neumDofs ) ) 

  elseif strcmp( modelProperties.analysisSettings.methodName, 'archLength' )
    
    %~ if norm(constantFext)>0 || ~(strcmp( userLoadsFilename , '')),
      %~ error('load case not implemented yet for Arc-Length method');
    %~ end
    
    %~ FextG  = computeFext( constantFext, variableFext, nextLoadFactor, userLoadsFilename ) ;
    
    %~ % incremental displacement
    %~ systemDeltauRHS = [ -(Fint(neumdofs)-FextG(neumdofs))  variableFext(neumdofs) ] ;

  elseif strcmp( modelProperties.analysisSettings.methodName, 'newmark' )

    FextG = computeFext( constantFext, variableFext, nextLoadFactor, userLoadsFilename ) ;
    
    rhat      =   Fint ( neumdofs ) ...
                + Fvis ( neumdofs ) ...
                + Fmas ( neumdofs ) ...
                - FextG( neumdofs ) ;
                
    systemDeltauRHS = -rhat ;

  elseif strcmp( modelProperties.analysisSettings.methodName, 'alphaHHT' )
      
    fs = assembler ( ...
      Conec, crossSecsParamsMat, coordsElemsMat, materialsParamsMat, KS, Ut, 1, Udott, ...
      Udotdott, nodalDispDamping, solutionMethod, elementsParamsMat ) ;
    
    Fintt = fs{1} ;  Fvist =  fs{2};  Fmast = fs{3} ;  

    FextG  = computeFext( constantFext, variableFext, nextLoadFactor, userLoadsFilename ) ;
    FextGt = computeFext( constantFext, variableFext, currLoadFactor, userLoadsFilename ) ;
                      
    rhat   =  ( 1 + alphaHHT ) * ( ...
                + Fint  ( neumdofs ) ...
                + Fvis  ( neumdofs ) ...
                - FextG ( neumdofs ) ...
              ) ...
              ...
              - alphaHHT * ( ...
                + Fintt ( neumdofs ) ...
                + Fvist ( neumdofs ) ...
                - FextGt( neumdofs ) ...
                ) ...
              ...
              + Fmas    ( neumdofs ) ;
                
    systemDeltauRHS = -rhat ;
    
  end
    
