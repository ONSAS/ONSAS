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
%%%%%%%%%%%%%% 	Matrices de masa con paso complejo   	%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MassMatrixComp,GyroMatrixComp] = MatrizMasaCompleja (xelem, Ue, Udote, Udotdote, params,Jrho )

	%-------------- Increment ----------------------%
	h  = 1e-10 ;

	%-------------- Mass Matrix Comp --------------%
	MassMatrixComp = zeros (12,12);
	VecIncrement = zeros (12,1)   ;

	for aux = 1:12
		VecIncrement(aux) = i*h;
		[FineMatrix] = elementFuerzaInercial (xelem, Ue, Udote, Udotdote+VecIncrement, params,Jrho);
		FineMatrix = 1/h * imag(FineMatrix);
		MassMatrixComp (:,aux) = FineMatrix;
		VecIncrement = zeros (12,1);
	end


	%-------------- Gyro Matrix Comp --------------%
    GyroMatrixComp = zeros (12,12);
	VecIncrement = zeros (12,1)   ;
	for aux = 1:12
		VecIncrement(aux) = i*h;
		[FineMatrix] = elementFuerzaInercial (xelem, Ue, Udote+VecIncrement, Udotdote, params,Jrho);
		FineMatrix = 1/h * imag(FineMatrix);
		GyroMatrixComp (:,aux) = FineMatrix;
		VecIncrement = zeros (12,1);
	end

end
