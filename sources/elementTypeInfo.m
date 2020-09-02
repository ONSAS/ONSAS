% Copyright (C) 2020, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
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


function [numNodes, dofsStep] = elementTypeInfo ( elemType )

switch elemType
case 1 % node
  numNodes = 1 ;
  dofsStep = 1 ;

case 2 % truss
  numNodes = 2 ;
  dofsStep = 2 ;

case 3
  numNodes = 2 ;
  dofsStep = 1 ;

case 4
  numNodes = 4 ;
  dofsStep = 2 ;

end
