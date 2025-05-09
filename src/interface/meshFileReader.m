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
function [nodesMat, conecMat] = meshFileReader(fileName)

  fileExtension = fileName((end - 2):end);

  if strcmp(fileExtension, 'msh')
    % md reads data from msh file
    [nodesMatinp, conecMatinp, physicalNames] = mshFormatReader(fileName);

    % md converts strings to integers indexes matrix
    matInds = zeros(length(physicalNames), 3);
    for i = 1:size(matInds, 1)
      for j = 1:3
        matInds(i, j) = str2num(physicalNames{i}(1 + (j - 1) * 3 + (1:2)));
      end
    end

    nodesMat = zeros(size(nodesMatinp, 1), 3);
    nodesMat(:, 1:3) = nodesMatinp(:, 1:3);

    conecMat = zeros(size(conecMatinp, 1), 3 + 4);
    conecMat(:, 3 + (1:4)) = conecMatinp(:, 1:4);

    % md adds MEB parameters to elements with params defined
    indsNZ = find(conecMatinp(:, 5));
    conecMat(indsNZ, 1:3) = matInds(conecMatinp(indsNZ, 5), :);

  else
    error('extension not implemented. Please report an issue.');
  end
