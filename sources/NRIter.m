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


% -------------------------------------
% newton raphson iteration
% -------------------------------------

KTred  = KT  ( neumdofs, neumdofs ) ;

Resred = FintGk(neumdofs) - FextG(neumdofs) ;

% incremental displacement
deltaured = KTred \ ( - Resred ) ;

normadeltau = norm( deltaured   )   ;
normaUk     = norm( Uk(neumdofs ) ) ;

Uk ( neumdofs ) = Uk(neumdofs ) + deltaured ;
