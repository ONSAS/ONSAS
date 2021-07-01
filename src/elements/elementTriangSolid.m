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

%md Function for computation of nodal forces and tangent stiffness matrix for 2D 3 nodes triangle element with linearElastic behavior. The dofs are assumed to be on x-y.
%md

function [ Finte, KTe, stress ] = elementTriangSolid( ...
  elemCoords, elemDisps, elemConstitutiveParams, paramOut, t, planeStateFlag )

  E  = elemConstitutiveParams(2) ;
  nu = elemConstitutiveParams(3) ;

  if planeStateFlag == 1
    C = E / (1-nu^2) * [ 1   nu  0           ; ...
                         nu  1   0           ; ...
                         0   0   (1-nu )/2 ] ;
                         
  elseif planeStateFlag == 2
    C = E * (1-nu) / ( (1+nu)*(1-2*nu) ) * ...
                       [ 1          nu/(1-nu)  0                   ; ...
                         nu/(1-nu)  1          0                   ; ...
                         0          0          (1-2*nu)/(2*(1-nu)) ] ;
  end


E
nu
  x = elemCoords(1:3:end)' ;  y = elemCoords(2:3:end)' ;
x
y
C

  A = 0.5 * det( [ ones(1,3) ; x' ; y' ] )  % element area

  assert( A>=0, 'Element with negative area, check connectivity.')

  B = 1 / (2*A) * [ y(2)-y(3)  0         0   y(3)-y(1)  0         0   y(1)-y(2)  0         0 ;
                    0          x(3)-x(2) 0   0          x(1)-x(3) 0   0          x(2)-x(1) 0 ;
                    x(3)-x(2)  y(2)-y(3) 0   x(1)-x(3)  y(3)-y(1) 0   x(2)-x(1)  y(1)-y(2) 0 ]

  KTe = B' * C * B * A * t

  strain = B * elemDisps ;
  stress = C * strain ;

  Finte = KTe * elemDisps ;
