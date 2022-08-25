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
 
% ==============================================================================
function [ fun ] = shapeFuns( x,y,z , derivOrder )

  if derivOrder == 0
    fun = zeros(4,1) ;
    fun(1) = x ;
    fun(2) = 1 - x - y - z ;
    fun(3) = z ;
    fun(4) = y ;

  elseif derivOrder == 1
    fun = zeros( 3, 4 ) ;
    fun(1,1) = 1 ;
    fun(1:3,2) = [ -1, -1, -1 ] ;
    fun(3,3) = 1.0 ;
    fun(2,4) = 1.0 ;
  end

% ==============================================================================
