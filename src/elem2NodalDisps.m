% function that constructs the assembled Fext vector for one BCtype 

function [ nonHomDiriVals, diriDofs, nonHomDiriDofs ] = elem2NodalDisps ( Conec, indBC, elements, boundaryConds, Nodes )
  
  elemsWithBC = find( Conec(:,3) == indBC ) ;
  
  diriDofs = [] ;
  nonHomDiriVals = [] ;
  nonHomDiriDofs = [] ;
  
  for elemInd = 1:length( elemsWithBC );

    elem            = elemsWithBC( elemInd ) ;
    
    nodesElem       = nonzeros( Conec (elem, 5:end ) ) ;
    
    elemType        = elements.elemType{ Conec(elem,2 )}  ;

    auxDofs = nodes2dofs( nodesElem, 6 ) ; auxDofs = auxDofs(:);
    
    % nodal loads
    % -----------
    if strcmp( elemType, 'node') ; % node
    
      impoDofs = boundaryConds.impoDispDofs{ indBC } ;
      impoVals = boundaryConds.impoDispVals{ indBC } ;
      locNonHomDofs = find( impoVals ) ;

      if ~isempty( locNonHomDofs)
        nonHomDiriDofs = [ nonHomDiriDofs; auxDofs(  locNonHomDofs) ];
        nonHomDiriVals = [ nonHomDiriVals; impoVals( locNonHomDofs) ];
      end

      diriDofs = [ diriDofs ; auxDofs(impoDofs) ] ;
      
    end %if elemTypes
    
  end % for elements
  
