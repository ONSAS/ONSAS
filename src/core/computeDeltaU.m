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
 
function [deltaured, nextLoadFactorVals ] = computeDeltaU( ...
  systemDeltauMatrix, systemDeltauRHS, dispIter, convDeltau, analysisSettings, nextLoadFactorVals, currDeltau, timeIndex, neumDofs )

ndofs = 6 ;
relDeltas = zeros(length(convDeltau)) ;
% relative disp
for i =1:( length(convDeltau)/ndofs -1 )
	dofsElem = ( (i-1)*ndofs+1:i*ndofs ) ;  
	redDeltaElem = convDeltau( (i*ndofs+1):(i+1)*ndofs ) - convDeltau( ((i-1)*ndofs+1):i*ndofs ) ;
	relDeltas(dofsElem) = abs(redDeltaElem) ;
end
relDeltas = relDeltas(neumDofs) ;

arcLengthNorm = zeros( size( convDeltau ) ) ;
arcLengthNorm(1:2:end) = 1 ;
arcLengthNorm = arcLengthNorm(neumDofs) ;

% keep reduced converged delta u
convDeltau = convDeltau( neumDofs ) ;

  if strcmp( analysisSettings.methodName, 'arcLength' )
		
    aux = systemDeltauMatrix \ systemDeltauRHS ;
    deltauast = aux(:,1) ;  deltaubar = aux(:,2) ;
    
    posVariableLoadBC = analysisSettings.posVariableLoadBC ;
		
		if length(analysisSettings.incremArcLen) > 1
			incremArcLen = analysisSettings.incremArcLen(timeIndex) ;
		else	
			incremArcLen = analysisSettings.incremArcLen ;
		end
		
		% =========================
		nodeControlBool = 0 ;
		cylindricalConstraintBool = 1 ;
		% =========================

  if dispIter == 1 % predictor solution
    if norm( convDeltau ) == 0
      deltalambda = analysisSettings.iniDeltaLamb ;
    else
      deltalambda = sign( convDeltau' * (arcLengthNorm .* deltaubar ) ) * incremArcLen / sqrt( deltaubar' * ( arcLengthNorm .* deltaubar ) ) ;
    end
	elseif nodeControlBool == 1 % nodal disp control
			
			%cMatrix = zeros(size(convDeltau,1),1) ;
			%dofs = nodes2dofs(node, gdlpernode) ;	
			%cMatrix(1:2:end) = 1 ;
			%cMatrix = cMatrix(neumDofs) ;
			%cMatrix = arcLengthNorm ;
			cMatrix = relDeltas ;
			deltalambda = (incremArcLen - cMatrix'*currDeltau - cMatrix'*deltauast ) / ( cMatrix'*deltaubar ) ;
  elseif cylindricalConstraintBool == 1  % cylindrical constraint equation
		% arcLengthNorm = relDeltas ;
    discriminant_not_accepted = true ;
    num_reductions = 0 ;

    while discriminant_not_accepted && num_reductions < 10

      ca =    deltaubar' * ( arcLengthNorm .* deltaubar) ;
      cb = 2*(currDeltau + deltauast)' * ( arcLengthNorm .* deltaubar ) ;
      cc =   (currDeltau + deltauast)' * ( arcLengthNorm .* (currDeltau + deltauast) ) - incremArcLen^2 ;
      disc = cb^2 - 4 * ca * cc ;

      if disc < 0
        cc
        disc
        num_reductions = num_reductions + 1 ;
        incremArcLen = incremArcLen * .5 ;
        warning( 'negative discriminant, reducing arc length time : %3i', num_reductions );
      else
        discriminant_not_accepted = false ;
      end

    end

    if disc < 0
      disc, error( 'negative discriminant');
    end

    sols = -cb/(2*ca) + sqrt(disc) / (2*ca)*[-1 +1]' ;

    % compute the scalar product
    vals = [ ( currDeltau + deltauast + deltaubar * sols(1) )' * ( arcLengthNorm .* currDeltau )   ;
              ( currDeltau + deltauast + deltaubar * sols(2) )' * ( arcLengthNorm .* currDeltau ) ] ;
    % choose lambda that maximices that scalar product
    deltalambda = sols( find( vals == max(vals) ) ) ;
  end
  
  nextLoadFactorVals( posVariableLoadBC )  = nextLoadFactorVals( posVariableLoadBC ) + deltalambda(1) ;

  deltaured = deltauast + deltalambda(1) * deltaubar ;

  else   % incremental displacement
    deltaured = systemDeltauMatrix \ systemDeltauRHS ;
  end
