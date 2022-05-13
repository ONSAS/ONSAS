function [ vtkNodes, vtkConec, vtkNodalDisps, vtkNormalForces ] ...
   = trussVtkData( Nodes, Conec, elemCrossSecParams, U )

  vtkNodes = [] ;
  vtkConec = [] ;
  vtkNodalDisps   = [] ;
  vtkNormalForces = [] ;

  nelem = size(Conec,1) ;

  counterNodes     = 0 ;

  for i=1:nelem
    nodesElem    = Conec(i,1:2) ;
    dofsElem     = nodes2dofs( nodesElem, 6 ) ;

    coordsElemNodes = reshape( Nodes( nodesElem(:), : )', 6, 1 ) ;

    % length of current element
    elemLength = norm( coordsElemNodes(4:6) - coordsElemNodes(1:3) ) ;

    % q section tengo
    [ iniNodes, midNodes, endNodes, sectPar ] = crossSectionVtkSolidConnec( elemCrossSecParams ) ;

    dispsElem          = U( dofsElem ) ;
    thetaLocIniSubElem = zeros(3,1)    ;
    thetaLocEndSubElem = zeros(3,1)    ;

    [ R0, Rr, locDisp ] = elementBeamRotData( coordsElemNodes, dispsElem ) ;

    dispLocIniSubElem  = zeros(3,1)    ;
    dispLocEndSubElem  = [ locDisp(1); 0; 0 ]    ;

    coordLocSubElem = [ 0; elemLength ] ;

    [ Nodesvtk, Conecvtk, Dispsvtk ] = vtkBeam2SolidConverter( coordsElemNodes, ...
       dispsElem, coordLocSubElem, dispLocIniSubElem, dispLocEndSubElem, thetaLocIniSubElem, thetaLocEndSubElem, sectPar, Rr, R0 ) ;

       Conecvtk( :, 2:end ) = Conecvtk(:,2:end)+counterNodes ;
       vtkNodes             = [ vtkNodes ;     Nodesvtk ]    ;
       vtkConec             = [ vtkConec ;     Conecvtk ]    ;
       vtkNodalDisps        = [ vtkNodalDisps; Dispsvtk ]    ;

       counterNodes = counterNodes + size( vtkNodes, 1 ) ;

  end % for plot element subdivision
