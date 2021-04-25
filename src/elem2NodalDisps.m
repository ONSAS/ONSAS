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

    auxDofs = nodes2dofs( nodesElem, 6 ) ; auxDofs = auxDofs(:)

    impoDofs = boundaryConds.imposDispDofs{ indBC } ;
    impoVals = boundaryConds.imposDispVals{ indBC } ;
    locNonHomDofs = find( impoVals ) ;
    
    % nodal loads
    % -----------
    if strcmp( elemType, 'node') ; % node

      if ~isempty( locNonHomDofs)
        nonHomDiriDofs = [ nonHomDiriDofs; auxDofs(  locNonHomDofs) ];
        nonHomDiriVals = [ nonHomDiriVals; impoVals( locNonHomDofs) ];
      end

      diriDofs = [ diriDofs ; auxDofs(impoDofs) ] ;

    elseif strcmp( elemType, 'triangle') ; % triangle
      
      for j=1:length(impoDofs)
        diriDofs = [ diriDofs ; auxDofs( impoDofs(j):6:end) ] ;
      end

    end %if elemTypes
    
  end % for elements
