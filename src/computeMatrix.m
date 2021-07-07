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
function systemDeltauMatrix = computeMatrix( Conec, elements, Nodes, materials, KS, analysisSettings, Uk, Udott, Udotdott, neumdofs, nodalDispDamping ) ;

  % computes static tangent matrix
  [ ~, ~, mats ] = assembler( Conec, elements, Nodes, materials, KS, Uk, Udott, Udotdott, analysisSettings, [0 0 1], nodalDispDamping ) ;

  KT      = mats{1} ;
  if strcmp( analysisSettings.methodName, 'newmark' ) || strcmp( analysisSettings.methodName, 'alphaHHT' )

    dampingMat = mats{2} ;
    massMat    = mats{3} ;

    global spitMatrices
    if spitMatrices == true
      KTred = KT(neumdofs,neumdofs);
      massMatred = massMat(neumdofs,neumdofs);
      save('-mat', 'output/matrices.mat', 'KTred','massMatred' );
      figure
      spy(full(KTred))
      figure
      spy(full(massMatred))
      stop
    end

  end

  if strcmp( analysisSettings.methodName, 'newtonRaphson' ) || strcmp( analysisSettings.methodName, 'arcLength' )

    systemDeltauMatrix = KT ( neumdofs, neumdofs ) ;

  elseif strcmp( analysisSettings.methodName, 'newmark' )

    alphaNM = analysisSettings.alphaNM ;
    deltaNM = analysisSettings.deltaNM ;
    deltaT  = analysisSettings.deltaT  ;

    systemDeltauMatrix =                                   KT(         neumdofs, neumdofs ) ...
                         + 1/( alphaNM * deltaT^2)       * massMat(    neumdofs, neumdofs ) ...
                         + deltaNM / ( alphaNM * deltaT) * dampingMat( neumdofs, neumdofs )  ;

  elseif strcmp( analysisSettings.methodName, 'alphaHHT' )

    alphaHHT = analysisSettings.alphaHHT ;
    deltaT   = analysisSettings.deltaT  ;

    deltaNM = (1 - 2 * alphaHHT ) / 2 ;
    alphaNM = (1 - alphaHHT ^ 2 ) / 4 ;

    systemDeltauMatrix = (1 + alphaHHT )                                 * KT         ( neumdofs, neumdofs ) ...
                       + (1 + alphaHHT ) * deltaNM / ( alphaNM*deltaT  ) * dampingMat ( neumdofs, neumdofs )  ...
                       +                         1 / ( alphaNM*deltaT^2) * massMat    ( neumdofs, neumdofs ) ;

  end
