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

massMat = sparse( ndofpnode*nnodes, ndofpnode*nnodes) ;

for elem = 1:nelems

  A  = secGeomProps(Conec(elem,6),1) ;
  rho = hyperElasParamsMat(Conec(elem,5),end) ;
  nodeselem = Conec(elem,1:2)' ;
  dofselem  = nodes2dofs( nodeselem , ndofpnode ) ;

  if Conec(elem,7)==1
    massMate = elementTruss3DMassMats(coordsElemsMat(elem,:), rho, A ) ;

  elseif Conec(elem,7)==2
    massMate = elementBeam3DMassMats( coordsElemsMat(elem,:), rho, A, deltamassMat );
  end
    % matrices assembly
  massMat (dofselem, dofselem ) = massMat (dofselem,dofselem) + massMate     ;
end
  % ------------------------------------
