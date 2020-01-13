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


% sets current tangent stiffness and geometric matrices
KL0    = sparse( 2*nnodes,2*nnodes) ;
KM     = sparse( 2*nnodes,2*nnodes) ;
KG     = sparse( 2*nnodes,2*nnodes) ;
FintGk = zeros( 2*nnodes,1) ;

for elem = 1:nelems

  % obtains nodes and dofs of element
  nodeselem = Conec(elem,1:2)' ;
  dofselem  = nodes2dofs( nodeselem , 2 ) ;
  
  % calculates important parameters of the element
  A  = As(Conec(elem,4));    le = largosini(elem);

  Xe = [ xelems( elem, 1) yelems( elem, 1) xelems( elem, 2) yelems( elem, 2) ]' ;
  Ue = Uk(dofselem) ;

  B1e = 1.0 / ( le^2 ) * Xe' * Ge' ;     B2e = 1.0 / ( le^2 ) * Ue' * Ge' ;
  
  KT1e  = dsigdepsk(elem) * A * le * ( B1e' * B1e ) ;
  KT2e  = dsigdepsk(elem) * A * le * ( B2e' * B1e + B1e' * B2e + B2e' * B2e ) ;
  Ksige = A * Stressk(elem) / le * Ge  ;
  Finte = (B1e+B2e)' * A * le * Stressk(elem) ;
  
  % matrices assembly
  KL0 (dofselem,dofselem) = KL0 (dofselem,dofselem) + KT1e        ;
  KM  (dofselem,dofselem) = KM  (dofselem,dofselem) + KT1e + KT2e ;
  KG  (dofselem,dofselem) = KG  (dofselem,dofselem) + Ksige       ;
  
  % internal loads vector assembly
  FintGk ( dofselem) = FintGk(dofselem) + Finte ;

end

KM  = KM  + KS ;
KL0 = KL0 + KS ;

FintGk = FintGk + KS * Uk ;

% boundary conditions are applied
KL0red = KL0 ( neumdofs, neumdofs );
KMred  = KM  ( neumdofs, neumdofs );
KGred  = KG  ( neumdofs, neumdofs );
KTred  = KMred + KGred ;
