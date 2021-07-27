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
%%%%%%%%%%% 			Matriz Kg del elemento	      	 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Kg_e] = elementMatrizKglobal (Dte,Ddote,Ddotdote,params,Jrho,xelem,BooleanComplexModes,NewmarkParams,numericalMethodParams)

	 beta   = NewmarkParams(1);
     alpha  = NewmarkParams(2);
	 DeltaT = numericalMethodParams(4);

	%loads Bt
	[Bt_e] = elementMatrizBt (Dte);

	%loads Mass and Gyroscopic Matrix
	if BooleanComplexModes ==1
		[MassMatrix_e,GyroMatrix_e] = MatrizMasaCompleja (xelem, Dte, Ddote, Ddotdote, params,Jrho );
	   else
		[Fine_e,MassMatrix_e,GyroMatrix_e] = elementFuerzaInercial(xelem, Dte, Ddote, Ddotdote, params,Jrho );
	end

	%Loads Ktang
	[ ~, KT_e] = elementBeam3DInternLoads(xelem, Dte, params ) ;

	%Compute Kg
	%~ MassMatrix_e
	%~ GyroMatrix_e
	Kg_e = + KT_e + 1/(beta*DeltaT^2) * MassMatrix_e * Bt_e + alpha/(beta*DeltaT) * GyroMatrix_e * Bt_e; %Ec 103 y 102

    %~ KG_din= + 1/(beta*DeltaT^2) * MassMatrix_e * Bt_e + alpha/(beta*DeltaT) * GyroMatrix_e * Bt_e
end
