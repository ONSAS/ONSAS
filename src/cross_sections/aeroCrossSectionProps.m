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


function [ chordVector, aeroCoefs ] = aeroCrossSectionProps ( elemCrossSecParams, chordVector, aeroCoefs)

    chordVecIsDefiend = ~isempty( chordVector ) ; 
    dragIsDefiend     = ~isempty( aeroCoefs{1}) ; 
    liftIsDefiend     = ~isempty( aeroCoefs{2}) ; 
    pitchIsDefiend    = ~isempty( aeroCoefs{3}) ; 


    if strcmp( elemCrossSecParams{1}, 'circle' ) || strcmp( elemCrossSecParams{1}, 'pipe' )

        if ~chordVecIsDefiend chordVector = [ 0 0 d ] end ;
        % if ~dragIsDefiend dragF = 'dragCoefCircular' ; end ;
        % if ~liftIsDefiend liftF = 'circularLift' ; end ;
        % if ~pitchIsDefiend pitchF = inline("0")  ; end ;
        if ~dragIsDefiend dragF = "dragCoefFunctionLA"  ; end ;
        if ~liftIsDefiend liftF = "liftCoefFunctionLA"  ; end ;
        if ~pitchIsDefiend pitchF = "liftCoefFunctionLA"; end ;
    end

    aeroCoefs = { dragF; liftF; pitchF } ;