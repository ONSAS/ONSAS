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
 
function [ Area, J, Iyy, Izz, Jrho ] = crossSectionProps ( elemCrossSecParams, rho )

if strcmp( elemCrossSecParams{1}, 'generic' ) %general section

    elemCrossSecParamsVec = elemCrossSecParams{2} ;
    Area = elemCrossSecParamsVec( 1 ) ;
    J    = elemCrossSecParamsVec( 2 ) ;
    Iyy  = elemCrossSecParamsVec( 3 ) ;
    Izz  = elemCrossSecParamsVec( 4 ) ;

    if length( elemCrossSecParamsVec ) > 5
        Jrho =  diag( elemCrossSecParamsVec( 5:7 ) ) ;
    else
        Jrho = rho * diag( [ J Iyy Izz ] ) ;
    end

elseif strcmp( elemCrossSecParams{1}, 'rectangle' )
    elemCrossSecParamsVec = elemCrossSecParams{2} ;
    Area = elemCrossSecParamsVec( 1 ) * elemCrossSecParamsVec( 2 )          ;
    Iyy  = elemCrossSecParamsVec( 1 ) * elemCrossSecParamsVec( 2 )^3 / 12.0 ;
    Izz  = elemCrossSecParamsVec( 2 ) * elemCrossSecParamsVec( 1 )^3 / 12.0 ;

    % torsional constant from table 10.1 from Roark's Formulas for Stress and Strain 7th ed.
    a = .5 * max( elemCrossSecParamsVec(1:2) ) ;
    b = .5 * min( elemCrossSecParamsVec(1:2) ) ;

    J = a * b^3 * ( 16/3 - 3.36 * b/a * ( 1 - b^4 / ( 12*a^4 ) ) ) ;

    Jrho = rho * diag( [ J Iyy Izz ] ) ;

elseif strcmp( elemCrossSecParams{1}, 'circle' )
    diameter = elemCrossSecParams{2} ;
    Area = pi*diameter^2/4             ;
    Iyy  = pi*diameter^4/64            ;
    Izz  = Iyy                         ;
    J    = Iyy + Izz                   ;
    Jrho = rho * diag( [ J Iyy Izz ] ) ;
    
elseif strcmp( elemCrossSecParams{1}, 'pipe' )
    elemCrossSecParamsVec = elemCrossSecParams{2} ;
    dext = elemCrossSecParamsVec(1)    ;
    dint = elemCrossSecParamsVec(2)    ;
    Area = pi*(dext^2 - dint^2)/4      ;
    Iyy  = pi *(dext^4 - dint^4) / 64  ;
    Izz  = Iyy                         ;
    J    = Iyy + Izz                   ;
    Jrho = rho * diag( [ J Iyy Izz ] ) ;
else
  error(' section type not implemented yet, please create an issue')
end
