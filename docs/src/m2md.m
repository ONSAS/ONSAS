% Copyright 2024, ONSAS Authors (see documentation)
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

function m2md(fileIn, fileOut, includeCodeBoolean, iniLine)

  fidIn  = fopen(fileIn, 'r');
  fidOut = fopen(fileOut, 'w');

  for i = 1:iniLine
    currentLine = fgetl(fidIn);
  end

  isInCodeBlock = ~(length(currentLine) >= 3 && strcmp(currentLine(1:3), '% md'));
  lineCount = 0;

  while ~feof(fidIn)

    lineCount = lineCount + 1;
    if lineCount ~= 1
      currentLine = fgetl(fidIn);
    end

    if length(currentLine) >= 7 && strcmp(currentLine((end - 6):end), '% hidden')
      % hidden line do not do anything

    elseif length(currentLine) >= 4 && strcmp(currentLine(1:4), '% md') % not code

      if isInCodeBlock % closes code block before writing comment
        if includeCodeBoolean
          fprintf(fidOut, '```\n');
        end
        isInCodeBlock = false;
      end
      fprintf(fidOut, '%s\n', currentLine(4:end));

    else
      if ~isInCodeBlock % open code block
        if includeCodeBoolean
          fprintf(fidOut, '```\n');
        end
        isInCodeBlock = true;
      end
      if includeCodeBoolean
        fprintf(fidOut, '%s\n', currentLine);
      end
    end
  end

  fclose(fidIn);
  fclose(fidOut);
