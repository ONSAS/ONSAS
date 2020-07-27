% Copyright (C) 2019, Jorge M. Perez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquin Viera, Mauricio Vanzulli  
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

function fileName = auxReadDXF(nomArch)

%~ cd '../input'
cd input
fid = fopen(nomArch);
cd '../sources'
fileName = 'auxTxt.txt' ;

auxFile = fopen(fileName, 'w+');

line = fgetl(fid);
s={};
while ischar(line)
	s=[s;line];
	if length(line) == 0
		line = 'lala';
	else
		line = fgetl(fid);
	end
end
fclose(fid);

rows = size(s,1);
for i=1:rows
	fprintf(auxFile, '%s\n', s{i});
end

cd ..
fclose(auxFile);
