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
% function for computation of the normal force and tanget matrix
% of 3D truss elements using engineering strain.


function Fther = elementTrussThermalForce( Xe, Ue, E, A, thermalExpansion, temperature )

  Xe    = Xe'     ;
  Xedef = Xe + Ue ;

  Bdif = [ -eye(3) eye(3) ] ;
  Ge   = Bdif' * Bdif       ;

  % initial/deformed lengths
  lini = sqrt( sum( ( Bdif * Xe    ).^2 ) ) ;
  ldef = sqrt( sum( ( Bdif * Xedef ).^2 ) ) ;

  % normalized reference and deformed co-rotational vector
  e1ref = Bdif * Xe    / lini ;
  e1def = Bdif * Xedef / ldef ;

  b1 = 1/(lini^2) * Xe' * Ge ;

  TTcl              = Bdif' * e1def ;

  Fther  =  thermalExpansion * temperature * E * A * TTcl ;

