% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% Joaquin Viera, J. Bruno Bazzano, Marcelo Forets, Jean-Marc Battini.
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

function B = swtichToONSASBase(A)

%---------------- Change of basis Le and Battini 2014 -> ONSAS matrix  -------------------
Pch = sparse (6,6);
Pch (1,1) = 1;
Pch (3,2) = 1;
Pch (5,3) = 1;
Pch (2,4) = 1;
Pch (4,5) = 1;
Pch (6,6) = 1;

P = [Pch sparse(6,6);sparse(6,6) Pch];

    if size(A,2)>1
		B = P*A*P';
	 else
		B = P*A;
	end
end