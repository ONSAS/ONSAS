% Copyright 2025, ONSAS Authors (see documentation)
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
function [deltaured, nextLoadFactorVals] = computeDeltaU( ...
                                                         systemDeltauMatrix, systemDeltauRHS, dispIter, convDeltau, analysisSettings, nextLoadFactorVals, currDeltau, timeIndex, neumDofs, args)

  global flag_predictor
  global arcLengthFlag % 1: cylindrical 2: jirasek

  if isempty(analysisSettings.ALdominantDOF)
    arcLengthFlag = 1;
  else
    arcLengthFlag = 2;
  end

  if strcmp(analysisSettings.methodName, 'arcLength')

    arcLengthNorm = args{1};
    incremArcLen = args{2};

    aux = systemDeltauMatrix \ systemDeltauRHS;

    deltauast = aux(:, 1);
    deltaubar = aux(:, 2);

    posVariableLoadBC = analysisSettings.posVariableLoadBC;

    if dispIter == 1 % predictor solution
      if norm(convDeltau) == 0
        deltalambda = analysisSettings.iniDeltaLamb;

      else
        deltalambda = sign(convDeltau' * (arcLengthNorm .* deltaubar)) * incremArcLen / sqrt(deltaubar' * (arcLengthNorm .* deltaubar));
        if flag_predictor == 1

          % Follow the sign of the predictor work increment (incremental work)
          Fext = systemDeltauMatrix * (arcLengthNorm .* deltaubar);
          sign(deltaubar' * Fext);
          deltalambda(1) = sign(deltaubar' * Fext) * incremArcLen / sqrt(deltaubar' * (arcLengthNorm .* deltaubar));

        elseif flag_predictor == 2

          % Follow the sign of the stiffness determinant
          detKT = det(systemDeltauMatrix);
          sign(detKT);
          deltalambda(1) = sign(detKT) * incremArcLen / sqrt(deltaubar' * (arcLengthNorm .* deltaubar));
        end

      end

    elseif arcLengthFlag == 2 % Jirasek approach
      cMatrix = zeros(size(convDeltau)); % Jirasek
      % cMatrix = zeros(max(neumDofs),1 ) ; % see Jirasek

      % Variables to be defined by user
      dominantDofs      = analysisSettings.ALdominantDOF(1);

      assert(length(dominantDofs) == 1);
      vecind = find(neumDofs == dominantDofs);
      dominantDofsInd = vecind(1);

      scalingProjection = analysisSettings.ALdominantDOF(2);

      % Projection matrix
      cMatrix(dominantDofsInd) = scalingProjection; % reduced projection matrix

      % cMatrix(dominantDofs) = scalingProjection ;
      % cMatrix = cMatrix(neumDofs) ; % reduced projection matrix

      deltalambda = (incremArcLen - cMatrix' * currDeltau - cMatrix' * deltauast) / (cMatrix' * deltaubar);

    elseif arcLengthFlag == 1  % Cylindrical constraint equation
      discriminant_not_accepted = true;
      num_reductions = 0;

      while discriminant_not_accepted && num_reductions < 10

        ca =    deltaubar' * (arcLengthNorm .* deltaubar);
        cb = 2 * (currDeltau + deltauast)' * (arcLengthNorm .* deltaubar);
        cc =   (currDeltau + deltauast)' * (arcLengthNorm .* (currDeltau + deltauast)) - incremArcLen^2;
        disc = cb^2 - 4 * ca * cc;

        if disc < 0
          num_reductions = num_reductions + 1;
          incremArcLen = incremArcLen * .5;
          warning('negative discriminant, reducing arc length time : %3i', num_reductions);
        else
          discriminant_not_accepted = false;
        end

      end

      if disc < 0
        error('negative discriminant');
      end

      sols = -cb / (2 * ca) + sqrt(disc) / (2 * ca) * [-1 +1]';

      % compute the scalar product
      vals = [(currDeltau + deltauast + deltaubar * sols(1))' * (arcLengthNorm .* currDeltau)
              (currDeltau + deltauast + deltaubar * sols(2))' * (arcLengthNorm .* currDeltau)];
      % choose lambda that maximizes that scalar product
      deltalambda = sols(find(vals == max(vals)));
    end

    nextLoadFactorVals(posVariableLoadBC)  = nextLoadFactorVals(posVariableLoadBC) + deltalambda(1);

    deltaured = deltauast + deltalambda(1) * deltaubar;

  else   % incremental displacement

    % systemDeltauRHS
    % full(systemDeltauMatrix)
    deltaured = systemDeltauMatrix \ systemDeltauRHS;
  end
