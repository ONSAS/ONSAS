function [ vtkNodes, vtkConec, vtkNodalDisps, vtkNormalForces ] ...
   = trussVtkData( Nodes, Conec, elemTypeGeom, U )

  vtkNodes = [] ;
  vtkConec = [] ;
  vtkNodalDisps   = [] ;
  vtkNormalForces = [] ;

  nelem = size(Conec,1) ;

  for i=1:nelem
    nodesElem    = Conec(i,1:2) ;
    dofsElem     = nodes2dofs( nodesElem, 6 ) ;
    coordSubElem = reshape( Nodes( nodesElem(:), : )', 6,1) ;

    % q section tengo
    [ iniNodes, midNodes, endNodes, sectPar ] = crossSectionVtkSolidConnec( elemTypeGeom ) ;

    nodalDisp = U( dofsElem(1:2:end) ) ;
    nodalRot  = zeros(6,1)             ;

    [ NodesCell, ConecCell ] = vtkBeam2SolidConverter ( coordSubElem, nodalDisp, nodalRot, sectPar ) ;

    % add number of previous nodes to the new nodes
    lastNode = size( vtkNodes, 1 ) ;
    ConecCell(:,2:end) = ConecCell(:,2:end) + lastNode ;

    vtkNodes = [ vtkNodes ; NodesCell ] ;
    vtkConec = [ vtkConec ; ConecCell ] ;
    vtkNodalDisps = [ vtkNodalDisps ; [ ones(4,1)*nodalDisp(1:3)'; ones(4,1)*nodalDisp(4:6)' ] ] ;

  end % for plot element subdivision
