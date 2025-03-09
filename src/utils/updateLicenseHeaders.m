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

% Recursive function for updating headers of source files in src and downwards ...

function updateLicenseHeaders(folder, root)

  lengthOfCurrentHeader = 17;

  if root == true && exist('newFileHeader.txt') ~= 2
    system(['head -n ' num2str(lengthOfCurrentHeader) ' updateLicenseHeaders.m > newFileHeader.txt']);
  end

  folder;
  files = dir(folder);

  for i = 1:length(files)
    i;
    if files(i).isdir
      if ~strcmp(files(i).name(1), '.')
        files(i).name;
        updateLicenseHeaders([folder '/' files(i).name], false);
      end
    else
      completeFilename = [folder '/' files(i).name];
      showHeaderAndReplace(completeFilename, lengthOfCurrentHeader);
    end
  end

  if root == true && exist('newFileHeader.txt') == 2
    system(['rm newFileHeader.txt']);
  end

function showHeaderAndReplace(filename, lengthOfCurrentHeader)

  system(['head -n ' num2str(lengthOfCurrentHeader) ' ' filename]);

  reply = input('    replace yes or no? (y/n):', 's');

  if strcmp(reply, 'y')
    system(['more currentFileHeader.txt > aux.txt']);
    system(['tail  -n +' num2str(lengthOfCurrentHeader + 1) ' ' filename ' >> aux.txt']);
    system(['mv aux.txt ' filename]);
    disp('file updated.');
  end
  pause;
