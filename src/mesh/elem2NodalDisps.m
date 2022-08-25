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
 
%md function that constructs the vectors of constrained degrees of freedom

function [ nonHomDiriVals, diriDofs, nonHomDiriDofs ] = elem2NodalDisps ( Conec, indBC, elemsWithBC, elements, impoDofs, impoVals, Nodes )

  % declare outputs
  nonHomDiriVals = [] ;
  diriDofs = [] ;
  nonHomDiriDofs = [] ;

  % find not null indexes of impoVals
  locNonHomDofs = find( impoVals )       ;

  %md loop in the elements to convert to nodal constraints
  for elemInd = 1:length( elemsWithBC );

    elem        = elemsWithBC( elemInd )             ;
    nodesElem   = nonzeros( Conec(elem, 5:end ) )    ;
    elemType    = elements( Conec(elem,2) ).elemType ;

    %md compute an auxiliar column vector with the global degrees of freedom of the nodes of the current element
    auxDofs = nodes2dofs( nodesElem, 6 ) ; auxDofs = auxDofs(:);

    %md nodal constraints
    if strcmp( elemType, 'node') ; % node
      if ~isempty( locNonHomDofs)
        nonHomDiriDofs = [ nonHomDiriDofs; auxDofs(  locNonHomDofs) ];
        nonHomDiriVals = [ nonHomDiriVals; impoVals( locNonHomDofs)' ];
      end
      diriDofs = [ diriDofs ; auxDofs(impoDofs) ] ;

    %md edge or triangle constraints
    elseif strcmp( elemType, 'triangle') || strcmp( elemType, 'edge')

      for j=1:length(impoDofs)
        diriDofs = [ diriDofs ; auxDofs( impoDofs(j):6:end) ] ;
      end

    end %if elemTypes

  end % for elements
