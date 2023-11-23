% Copyright 2023, ONSAS Authors (see documentation)
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
function [ chordVector, aeroCoefs ] = aeroCrossSectionProps ( elemCrossSecParams, chordVector, aeroCoefs)

    chordVecIsDefined = ~isempty( chordVector ) ; 
    dragIsDefined     = ~isempty( aeroCoefs{1}) ; 
    liftIsDefined     = ~isempty( aeroCoefs{2}) ; 
    pitchIsDefined    = ~isempty( aeroCoefs{3}) ; 

    if strcmp( elemCrossSecParams{1}, 'circle' ) || strcmp( elemCrossSecParams{1}, 'pipe' )

        d_ext = elemCrossSecParams{2} ;
        if ~chordVecIsDefined chordVector = [ 0 0 d_ext ] ; end ;
        if ~dragIsDefined dragF = 'innerDragCoefCircular'; else dragF = aeroCoefs{1};  end ;
        anonymus_null = @(beta,Re) 0 ;
        if ~liftIsDefined liftF = anonymus_null ; else liftF = aeroCoefs{2}; end ;
        if ~pitchIsDefined pitchF = anonymus_null  ; else pitchF = aeroCoefs{3};  end ;

    end

    aeroCoefs = { dragF; liftF; pitchF } ;