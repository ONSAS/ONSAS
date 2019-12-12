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

%function for computation of the normal force of 3D truss elements

function [Finte, stress, dstressdeps, strain] = elementTruss3DInternLoads( ...
  Xe, Ue, hyperelasparams, A )

  localAxisRef = Xe( [(1:2:5)+6]) - Xe( [(1:2:5)]);
  lini = sqrt( sum( localAxisRef.^2 ) ) ;
  
  Xedef = Xe(:) + Ue(:) ;
  localAxisDef = Xedef( [(1:2:5)+6]) - Xedef( [(1:2:5)]);
  ldef = sqrt( sum( localAxisDef.^2 ) ) ;

  e1def = localAxisDef / ldef ;

  strain = ( ldef^2 - lini^2 ) / ( lini * (lini + ldef) ) ;

  [ stress, dstressdeps  ] = hyperElasModels ( strain, hyperelasparams ) ;

  TTcl              = zeros(12,1) ;
  TTcl(1:2:5)       = -e1def ;
  TTcl([(1:2:5)+6]) =  e1def ;

  Finte = stress * A * TTcl  ;
