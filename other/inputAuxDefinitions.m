% Copyright (C) 2019, Jorge M. Perez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquin Viera, Mauricio Vanzulli  
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


% ---------------------------------------------------


%~ releasesDofs = [];
%~ for i=1:nelems
  %~ if Conec(i,7)==1
    %~ releasesDofs = [ releasesDofs; nodes2dofs( Conec(i,1:2),ndofpnode)(2:2:end) ];
    %~ Releases = [ Releases; i ones(1,4) ] ;
  %~ end
%~ end

%~ releasesDofs = unique( releasesDofs);


% -------------------- indexes computation --------------------
indexesElems = zeros(nelems,1) ;
ntruss = 0 ;
nbeam = 0 ;
ntet = 0 ;
nplate = 0 ;
trussElem = [] ;
beamElem = [] ;
beamNodes = [] ;
tetElem = [] ;
plateElem = [] ;

for i = 1:nelems
  if Conec(i,7) == 1
    nodeselem = Conec(i,1:2) ;
    ntruss = ntruss + 1 ;
    indexesElems(i) = nbeam+1 ;
    trussElem = [ trussElem ; i ] ;
  elseif Conec(i,7) == 2
    nodeselem = Conec(i,1:2) ;
    nbeam = nbeam + 1 ;
    indexesElems(i) = nbeam ;
    beamElem = [ beamElem ; i ] ;
    beamNodes = [ beamNodes ; nodeselem' ] ;
  elseif Conec(i,7) == 3
    ntet = ntet + 1 ;
    indexesElems(i) = ntet ;
    tetElem = [ tetElem ; i ] ;
  elseif Conec(i,7) == 4
    nplate = nplate + 1 ;
    indexesElems(i) = nplate ;
    plateElem  = [ plateElem ; i ] ;
  end
end

% ------------------------------------------------------------
