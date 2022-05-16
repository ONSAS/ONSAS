% Copyright 2022, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
% Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini. 

% This file is part of ONSAS 

% ONSAS is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or 
% (at your option) any later version. 

% ONSAS is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.
 
function [deltaured, nextLoadFactorVals ] = computeDeltaU ( systemDeltauMatrix, systemDeltauRHS, dispIter, redConvDeltau, analysisSettings, nextLoadFactorVals, currDeltau  )

  convDeltau = redConvDeltau ;

  if strcmp( analysisSettings.methodName, 'arcLength' )

    aux = systemDeltauMatrix \ systemDeltauRHS ;

    incremArcLen = analysisSettings.incremArcLen ;

    deltauast = aux(:,1) ;  deltaubar = aux(:,2) ;

    posVariableLoadBC = analysisSettings.posVariableLoadBC ;

    if dispIter == 1 % predictor solution
      if norm(convDeltau)==0
        deltalambda = analysisSettings.iniDeltaLamb ;
      else
        deltalambda = sign( convDeltau' * deltaubar ) * incremArcLen / sqrt( deltaubar' * deltaubar ) ;
      end

    else % cylindrical constraint equation
      ca =    deltaubar' * deltaubar ;
      cb = 2*(currDeltau + deltauast)' * deltaubar ;
      cc =   (currDeltau + deltauast)' * (currDeltau + deltauast) - incremArcLen^2 ;
      disc = cb^2 - 4 * ca * cc ;
      if disc < 0
        disc, error( 'negative discriminant');
      end
      sols = -cb/(2*ca) + sqrt(disc) / (2*ca)*[-1 +1]' ;

      % compute the scalar product
      vals = [ ( currDeltau + deltauast + deltaubar * sols(1) )' * currDeltau   ;
               ( currDeltau + deltauast + deltaubar * sols(2) )' * currDeltau ] ;
      % choose lambda that maximices that scalar product
      deltalambda = sols( find( vals == max(vals) ) ) ;
    end

    nextLoadFactorVals( posVariableLoadBC )  = nextLoadFactorVals( posVariableLoadBC ) + deltalambda(1) ;

    deltaured = deltauast + deltalambda(1) * deltaubar ;

  else   % incremental displacement
    deltaured = systemDeltauMatrix \ systemDeltauRHS ;

  end
