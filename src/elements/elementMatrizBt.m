% Copyright (C) 2021, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% 			Matriz Bt del elemento	      	 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Bt_e] = elementMatrizBt (Dte)

Bt 		 = zeros(size(Dte,1),size(Dte,1));

%Angulos en el paso k
Wtk = Dte(2:2:end);
%Matriz Bt
for aux= 1:6:size(Dte,1)
	Bt(aux:aux+2,aux:aux+2) = eye(3,3);

	%calculo de la Ts
	indexANG1 = aux/2 + 1/2 ;
	indexANG2	 = indexANG1 +2;
	Wti = Wtk(indexANG1:indexANG2);
	inversaTs=invTs(Wti);
	Bt(aux+3:aux+5,aux+3:aux+5) = inversaTs';
end

Bt_e = Cambio_Base(Bt);
