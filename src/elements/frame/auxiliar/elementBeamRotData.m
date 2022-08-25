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
 
function  [ R0, Rr, locDisp ] = elementBeamRotData( xs, Ue ) ;

  assert( iscolumn( xs ), 'the coordinates must be provided in a column vector.' )

  % global thetas
  tg1 = Ue( 2:2:6  ) ;
  tg2 = Ue( 8:2:12 ) ;

  % rotation matrices
  Rg1 = expon( tg1 ) ;
  Rg2 = expon( tg2 ) ;

  x21 = xs(4:6)    - xs(1:3)   ;
  d21 = Ue(7:2:11) - Ue(1:2:5) ;

  refLength = norm( x21       ) ;
  defLength = norm( x21 + d21 ) ;

  % rotation matrix to reference configuration
  R0 = beamRefConfRotMat( x21 ) ;

  % deformed x axis
  e1 = ( x21 + d21 ) / defLength   ;

  q1 = Rg1 * R0 * [0 1 0]'  ;
  q2 = Rg2 * R0 * [0 1 0]'  ;
  q  = ( q1 + q2 ) / 2      ;

  % deformed z local axis
  e3 = cross ( e1, q ) ;
  e3 = e3 / norm( e3 ) ; % normalization

  % deformed y local axis
  e2 = cross ( e3, e1 );

  % rotation matrix
  Rr = [ e1 e2 e3 ] ;
  % -------------------

  % --- local displacements ---
  % axial displacement
  ul  = defLength - refLength ;

  % local rotations
  % Rr * Re1 * u = Rg1 * R0 * u
  Re1 = Rr' * Rg1 * R0 ;
  Re2 = Rr' * Rg2 * R0 ;

  tl1 = logar( Re1 ) ;
  tl2 = logar( Re2 ) ;

  locDisp = [ ul tl1' tl2' ] ;
