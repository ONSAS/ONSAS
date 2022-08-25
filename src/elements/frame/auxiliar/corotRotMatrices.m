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

function [R0, Rr, Rg1, Rg2, Rroof1, Rroof2] = corotRotMatrices(Ue, elemCoords ) 

  
  % element coordinates
  xs = elemCoords(:) ;
  

  % extract displacements and write them using Battini's nomenclature 
  permutIndxs = [ 1:2:5 2:2:6 ([1:2:5] + 6) ([2:2:6] + 6) ] ;
  dg          = Ue( permutIndxs ) ;
    

  % -------- global rotation matrices Rgs ------------  
  % global thetas
  tg1 = dg(  4:6  ) ;
  tg2 = dg( 10:12 ) ;

  % rotation matrices
  Rg1 = expon( tg1 ) ;
  Rg2 = expon( tg2 ) ;
  % --------------------------------------------------


  % - coordinates and reference rotation Matrix R0 ---
  [x21, d21, l, l0] = corotLenCoords(xs ,dg) ;
  % rotation matrix to reference configuration
  R0 = beamRefConfRotMat( x21 ) ;
  % ---------------------------------------------------


  % -------------- rigid rotation matrix Rr ------------
  % deformed x axis
  e1 = ( x21 + d21 ) / l   ;
  % auxiliary vectors q
  q1 = Rg1 * R0 * [0 1 0]' ;
  q2 = Rg2 * R0 * [0 1 0]' ;
  q  = ( q1 + q2 ) / 2     ;

  % deformed z local axis
  e3 = cross(e1, q)    ;
  e3 = e3 / norm( e3 ) ; % normalization

  % deformed y local axis
  e2 = cross (e3, e1) ;

  % rotation matrix
  Rr = [ e1 e2 e3 ] ;
  % ---------------------------------------------------

  % ---------- local rotation matrix Rroof ------------
  % Rr * Re1 * u = Rg1 * R0 * u
  Rroof1 = Rr' * Rg1 * R0 ;
  Rroof2 = Rr' * Rg2 * R0 ;
  % ---------------------------------------------------

