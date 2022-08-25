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
 
function [ vtkNodes, vtkConec, vtkNodalDisps, vtkNormalForces ] ...
   = frameVtkData( Nodes, Conec, elemCrossSecParams, U )

  vtkNodes        = [] ;
  vtkConec        = [] ;
  vtkNodalDisps   = [] ;
  vtkNormalForces = [] ;

  nPlotSubElements = 10 ; % number of plot subsegments
  counterNodes     = 0 ;
  nelem            = size(Conec,1) ;

  for i=1:nelem

    % nodes and degrees of freedom of current element
    nodesElem  = Conec(i,1:2)               ;
    dofsElem   = nodes2dofs( nodesElem, 6 ) ;

    % column vector with the coordinates of the nodes of the element
    coordsElemNodes = reshape( Nodes( nodesElem(:), : )', 6, 1 ) ;

    % computes connectivity of vtk element associated with frame cross-section
    [ iniNodes, midNodes, endNodes, sectPar ] = crossSectionVtkSolidConnec( elemCrossSecParams ) ;

    % length of current element
    elemLength = norm( coordsElemNodes(4:6) - coordsElemNodes(1:3) ) ;

    % column vector with displacements of the dofs of the current element
    dispsElem  = U( dofsElem ) ;

    % column vector with discretization in subelements, using local coordinates (r1,r2,r3)
    xsloc  = linspace( 0, elemLength, nPlotSubElements+1 )' ;

    % interpolation functions evaluation matrix with columns
    %    [  N1(lin1)                          N2(lin2)            N3(v1) N4(theta1) N5(v2) N6(theta2) ]
    interFuncLinear = [ ( elemLength - xsloc )/elemLength  xsloc/elemLength  ] ;
    interFuncCubic  = bendingInterFuns( xsloc, elemLength, 0 ) ;
    interFuncQuad   = bendingInterFuns( xsloc, elemLength, 1 ) ;

    [ R0, Rr, locDisp ] = elementBeamRotData( coordsElemNodes, dispsElem ) ;

    ul              = locDisp( 1   ) ;
    thetaLocIniElem = locDisp( 2:4 ) ;
    thetaLocEndElem = locDisp( 5:7 ) ;

    % interpolation of displacements
    valsLocDispXSubElements = interFuncLinear * [ 0;  ul ] ;
    valsLocDispYSubElements = interFuncCubic  * [ 0;  thetaLocIniElem(3); 0;  thetaLocEndElem(3) ] ;
    valsLocDispZSubElements = interFuncCubic   * [ 0; -thetaLocIniElem(2); 0; -thetaLocEndElem(2) ] ;

    % interpolation of angles
    valsLocThetaXSubElements = interFuncLinear * [    thetaLocIniElem(1);    thetaLocEndElem(1) ] ;
    valsLocThetaYSubElements = interFuncQuad   * [ 0; thetaLocIniElem(2); 0; thetaLocEndElem(2) ] ;
    valsLocThetaZSubElements = interFuncQuad   * [ 0; thetaLocIniElem(3); 0; thetaLocEndElem(3) ] ;

    for j = 1:nPlotSubElements,

      dispLocIniSubElem = [ valsLocDispXSubElements( j   ) ; ...
                            valsLocDispYSubElements( j   ) ; ...
                            valsLocDispZSubElements( j   ) ] ;

      dispLocEndSubElem = [ valsLocDispXSubElements( j+1 ) ; ...
                            valsLocDispYSubElements( j+1 ) ; ...
                            valsLocDispZSubElements( j+1 ) ] ;

      % 2-compute local rotations for the initial and final sections
      thetaLocIniSubElem = [ valsLocThetaXSubElements( j   ) ; ...
                             valsLocThetaYSubElements( j   ) ; ...
                             valsLocThetaZSubElements( j   ) ] ;
      %
      thetaLocEndSubElem = [ valsLocThetaXSubElements( j+1 ) ; ...
                             valsLocThetaYSubElements( j+1 ) ; ...
                             valsLocThetaZSubElements( j+1 ) ] ;

      coordLocSubElem = [ xsloc(j); xsloc(j+1) ] ;

      [ Nodesvtk, Conecvtk, Dispsvtk ] = vtkBeam2SolidConverter( coordsElemNodes, ...
         dispsElem, coordLocSubElem, dispLocIniSubElem, dispLocEndSubElem, thetaLocIniSubElem, thetaLocEndSubElem, sectPar, Rr, R0 ) ;

      Conecvtk( :, 2:end ) = Conecvtk(:,2:end)+counterNodes ;
      vtkNodes             = [ vtkNodes ;     Nodesvtk ] ;
      vtkConec             = [ vtkConec ;     Conecvtk ] ;
      vtkNodalDisps        = [ vtkNodalDisps; Dispsvtk ] ;

      counterNodes = counterNodes + (size(Conecvtk,2)-1) ;

    end % for plot points

  end % for elements


  %
  %
  %       if j==1 % compute the first face of the first subelement
  %         % compute the transversal displacements ux,uy,uz for the baricenter of the
  %         % initial cross section of the current sub element
  %         % interpolates using navier-Bernoulli
  %         veculIni = [ 0 ; ...
  %                  [ Ns(j, 2+2)  Ns(j, 2+4) ]*[tl1(3); tl2(3)] ; ...
  %                 -[ Ns(j, 2+2)  Ns(j, 2+4) ]*[tl1(2); tl2(2)]  ]
  % Rr
  % absolutedisp = Rr        * veculIni
  % %stop
  %         % compute nodal displacemets, interpolates using linear functions
  %         nodalDispIni   = Ns(j  ,1) * ( dispsElem(1:2:6 ) ) ...
  %                        + Ns(j  ,2) * ( dispsElem(7:2:12) ) ...
  %                        + Rr        * veculIni
  %         coordPlotSubElement(1:3) = Ns(j  ,1) * ( coordSubElem(1:3 ) ) ...
  %                                  + Ns(j  ,2) * ( coordSubElem(4:6 ) )
  %
  %       else
  %         veculIni = veculEnd ;
  %         nodalDispIni = nodalDispEnd ;
  %         coordPlotSubElement(1:3) = coordPlotSubElement(4:6) ;
  %       end
  %
  %       veculEnd = [ 0 ; ...
  %                [ Ns(j+1, 2+2)  Ns(j+1, 2+4) ]*[tl1(3); tl2(3)] ; ...
  %               -[ Ns(j+1, 2+2)  Ns(j+1, 2+4) ]*[tl1(2); tl2(2)]  ] ;
  %
  %       nodalDispEnd   = Ns(j+1,1) * ( dispsElem(1:2:6 ) ) ...
  %                      + Ns(j+1,2) * ( dispsElem(7:2:12) ) ...
  %                      + Rr        * veculEnd ;
  %       coordPlotSubElement(4:6) = Ns(j+1,1) * ( coordSubElem(1:3 ) ) ...
  %                                + Ns(j+1,2) * ( coordSubElem(4:6 ) ) ;
  %
  %       nodalDisp = [ nodalDispIni ; nodalDispEnd ] ;
