%~ Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli  

%~ This file is part of ONSAS.

%~ ONSAS is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation, either version 3 of the License, or
%~ (at your option) any later version.

%~ ONSAS is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%~ GNU General Public License for more details.

%~ You should have received a copy of the GNU General Public License
%~ along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

% --------------------------------------------------------------------------------------------------

% ==============================================================================
function BMat = BMats ( deriv )

  BMat = zeros(6,12) ;
  
  for k = 1:4

    for i = 1:3
      BMat ( i , (k-1)*3 + i  ) = deriv(i,k) ;
    end

    BMat ( 4 , (k-1)*3 + 2   ) = deriv(3,k) ;
    BMat ( 4 , (k-1)*3 + 3   ) = deriv(2,k) ;

    BMat ( 5 , (k-1)*3 + 1   ) = deriv(3,k) ;
    BMat ( 5 , (k-1)*3 + 3   ) = deriv(1,k) ;

    BMat ( 6 , (k-1)*3 + 1   ) = deriv(2,k) ;
    BMat ( 6 , (k-1)*3 + 2   ) = deriv(1,k) ;
      
  end
end
% ==============================================================================
