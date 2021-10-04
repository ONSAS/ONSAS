function [ vtkNodes, vtkConec, vtkNodalDisps, vtkNormalForces ] ...
   = frameVtkData( Nodes, Conec, elemTypeGeom, U )

  vtkNodes = [] ;
  vtkConec = [] ;
  vtkNodalDisps   = [] ;
  vtkNormalForces = [] ;

  nPlotSubElements = 10 ; % number of plot subsegments
  counterNodes     = 0 ;
  nelem = size(Conec,1) ;

  for i=1:nelem

    nodesElem    = Conec(i,1:2)               ;
    dofsElem     = nodes2dofs( nodesElem, 6 ) ;
    coordSubElem = reshape( Nodes( nodesElem(:), : )', 6, 1 ) ;

    l = norm( coordSubElem(4:6) - coordSubElem(1:3) ) ; % length

    [ iniNodes, midNodes, endNodes, sectPar ] = crossSectionVtkSolidConnec( elemTypeGeom ) ;
    dispsElem  = U( dofsElem ) ;

    [ ~, ~, ~, rotData ] = elementBeamForces( coordSubElem, ones(5, 1)', [1 1 1], dispsElem, [], [] ,0 ) ;

    locDisp = rotData{1} ;
    Rr      = rotData{2} ;

    % set rotations in co-rot system
    ul  = locDisp(1)   ;   tl1 = locDisp(2:4) ;  tl2 = locDisp(5:7) ;

    % loc axial coords to evaluate magnitudes
    xsloc  = linspace( 0 , l, nPlotSubElements+1 )' ;

    % conecSubElem = zeros( nPlotPoints, 2 ) ;

    % interpolation functions evaluation  matrix with columns
    %    [  N1(lin1)        N2(lin2)   N3(v1) N4(theta1) N5(v2) N6(theta2) ]
    Ns = [ ( l - xsloc )/l  xsloc/l    bendingInterFuns( xsloc, l, 0 )     ] ;

    coordPlotSubElement = zeros(6,1) ;

    for j=1:nPlotSubElements,
      %fprintf('subelement: %3i \n--------- \n',j)

      if j==1
        veculIni = [ 0 ; ...
                 [ Ns(j, 2+2)  Ns(j, 2+4) ]*[tl1(3); tl2(3)] ; ...
                -[ Ns(j, 2+2)  Ns(j, 2+4) ]*[tl1(2); tl2(2)]  ] ;

        nodalDispIni   = Ns(j  ,1) * ( dispsElem(1:2:6 ) ) ...
                       + Ns(j  ,2) * ( dispsElem(7:2:12) ) ...
                       + Rr        * veculIni ;
        coordPlotSubElement(1:3) = Ns(j  ,1) * ( coordSubElem(1:3 ) ) ...
                                 + Ns(j  ,2) * ( coordSubElem(4:6 ) ) ;

      else
        veculIni = veculEnd ;
        nodalDispIni = nodalDispEnd ;
        coordPlotSubElement(1:3) = coordPlotSubElement(4:6) ;
      end

      veculEnd = [ 0 ; ...
               [ Ns(j+1, 2+2)  Ns(j+1, 2+4) ]*[tl1(3); tl2(3)] ; ...
              -[ Ns(j+1, 2+2)  Ns(j+1, 2+4) ]*[tl1(2); tl2(2)]  ] ;

      nodalDispEnd   = Ns(j+1,1) * ( dispsElem(1:2:6 ) ) ...
                     + Ns(j+1,2) * ( dispsElem(7:2:12) ) ...
                     + Rr        * veculEnd ;
      coordPlotSubElement(4:6) = Ns(j+1,1) * ( coordSubElem(1:3 ) ) ...
                               + Ns(j+1,2) * ( coordSubElem(4:6 ) ) ;

      nodalDisp = [ nodalDispIni ; nodalDispEnd ] ;

      [ Nodesvtk, Conecvtk ] = vtkBeam2SolidConverter( coordPlotSubElement, nodalDisp, locDisp(2:7), sectPar ) ;

      Conecvtk(:,2:end) = Conecvtk(:,2:end)+counterNodes ;
      vtkNodes = [ vtkNodes ; Nodesvtk ] ;
      vtkConec = [ vtkConec ; Conecvtk ] ;
      vtkNodalDisps = [ vtkNodalDisps; ones(4, 1)*(nodalDispIni') ; ones(4, 1)*(nodalDispEnd') ] ;

      counterNodes = counterNodes + 8 ;

    end % for plot points

  end % for elements
