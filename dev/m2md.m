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

  % Skip initial lines if specified
  for i = 1:iniLine
    currentLine = fgetl(fidIn);
  end

  isInCodeBlock = false;
  lineCount = 0;

  while ~feof(fidIn)
    lineCount = lineCount + 1;
    if lineCount ~= 1
      currentLine = fgetl(fidIn);
    end

    % Skip empty lines
    if isempty(currentLine)
      continue
    end

    % Handle markdown comments (lines starting with % md)
    if length(currentLine) >= 4 && strcmp(currentLine(1:4), '% md')
      if isInCodeBlock && includeCodeBoolean
        fprintf(fidOut, '```\n'); % Close code block before markdown
        isInCodeBlock = false;
      end
      % Write markdown content without the '% md' prefix
      fprintf(fidOut, '%s\n', strtrim(currentLine(5:end)));

      % Handle hidden comments
    elseif length(currentLine) >= 8 && strcmp(currentLine((end - 7):end), '% hidden')
      % Skip hidden lines
      continue

      % Handle regular code
    else
      % Skip pure comment lines that aren't markdown
      if length(currentLine) >= 1 && currentLine(1) == '%'
        continue
      end

      if includeCodeBoolean
        if ~isInCodeBlock
          fprintf(fidOut, '```\n'); % Open new code block
          isInCodeBlock = true;
        end
        fprintf(fidOut, '%s\n', currentLine);
      end
    end
  end

  % Close final code block if needed
  if isInCodeBlock && includeCodeBoolean
    fprintf(fidOut, '```\n');
  end

  fclose(fidIn);
  fclose(fidOut);
