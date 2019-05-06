%~ Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli  

%~ This file is part of ONSAS.

%~ ONSAS is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation, either version 3 of the License, or
%~ (at your option) any later version.

%~ ONSAS is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%~ GNU General Public License for more details.

%~ You should have received a copy of the GNU General Public License
%~ along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.


% script for vtk output data generation

xs = Nodes(:,1) ;
ys = Nodes(:,2) ;
zs = Nodes(:,3) ;


function [sk] = skew(x);
sk = [   0 -x(3)  x(2) ;
      x(3)    0  -x(1) ;
     -x(2)  x(1)     0 ];
end

function R = expon(t);
  I = eye(3,3);  al = norm(t) ;
  if al == 0
    R = I ;
  else
    Rsk = skew(t);
    R = I+sin(al)/al*Rsk+2*(sin(al/2)/al)^2*Rsk^2;
  end
end

% ------------------------------------------------------------------------------
function [Conecvtk, Nodesvtk, vtkDispMat, vtkNormalForceMat ] = vtkConecNodes ( Nodes, Conec, secc, UG, nonLinearAnalysisBoolean, dynamicAnalysisBoolean, coordsElem, normalForce )
% ------------------------------------------------------------------------------

  nelems = size(Conec,1) ; nnodes = size(Nodes,1) ;
  ndofpnode = 6 ;
  %
  dispsElemsMat = zeros(nelems,2*ndofpnode) ;
  Local2GlobalMatsDef = cell(nelems,1) ;

	NodesAux                = [] ;
	NodesDefAux             = [] ;
  RotsAux                 = [] ;
	ConecAux                = [] ;
  vtkNormalForceMat       = [] ; 

  for i = 1:nelems
		elemMat = Conec(i,5) ;
    elemSec = Conec(i,6) ;
    elemType = Conec(i,7) ;
    
    if elemType == 1 || elemType == 2 && secc(1)>0
    
      nodeselem = Conec(i,1:2)' ;
      dofselem  = nodes2dofs( nodeselem , ndofpnode ) ;
      dispsElemsMat( i, : ) = UG(dofselem)' ;
      
      if nonLinearAnalysisBoolean == 1 || dynamicAnalysisBoolean == 1
      
        [xdef, ydef, zdef, conecElem] = outputFrameElementPlot( coordsElem(i,:)', dispsElemsMat(i,:)', elemType ) ;
        
      elseif nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
    
        [~, locglos] = beamParameters( Nodes(nodeselem,:) ) ;
        Local2GlobalMatsDef{i} = locglos ;
        [xdef, ydef, zdef, titax, titay, titaz, conecElem] = outputFrameElementPlotLin( coordsElem(i,:)', dispsElemsMat(i,:)', elemType, locglos ) ; 
        
      end

      % Element division 
      if elemType == 2
        ndivselem = conecElem(end,2) ;
      else
        ndivselem = 2 ;
      end
      % Undeformed nodes
      xloc = linspace( Nodes(nodeselem(1),1), Nodes(nodeselem(2),1), ndivselem)' ;  
      yloc = linspace( Nodes(nodeselem(1),2), Nodes(nodeselem(2),2), ndivselem)' ;  
      zloc = linspace( Nodes(nodeselem(1),3), Nodes(nodeselem(2),3), ndivselem)' ;
      % Undeformed nodes matrix
      NodesAux 		= [ NodesAux			; xloc yloc zloc ] ;   
      % Deformed nodes matrix
      NodesDefAux = [ NodesDefAux		; xdef ydef zdef ] ;
      % Rots aux
      RotsAux     = [ RotsAux       ; titax titay titaz ] ;   
      % New connectivity
      if elemType == 2
        ConecAux = [ ConecAux ; ((i-1)*ndivselem+conecElem) zeros(size(conecElem,1),2) ...
                                ones(size(conecElem,1),1)*elemMat ones(size(conecElem,1),1)*elemSec ones(size(conecElem,1),1)*elemType ] ; 

      else 
        ConecAux = [ ConecAux ; 2*i-1 2*i zeros(size(conecElem,1),2) ...
                                ones(size(conecElem,1),1)*elemMat ones(size(conecElem,1),1)*elemSec ones(size(conecElem,1),1)*elemType ] ;  
      end

      vtkNormalForceMat = [ vtkNormalForceMat ; normalForce(i)*ones(ndivselem-1,1) ] ;
    end
  end
  
  if size(ConecAux,1) ~= 0  
    % New UG
    UG = zeros( size(NodesAux,1)*ndofpnode, 1 ) ;
    % Nodes dispMat for VTK
    dispMat = NodesDefAux - NodesAux ;
    % Replaces nodes x, y, z disps 
    UG(1:2:end) = reshape( dispMat', size(NodesAux,1)*3, 1) ;
    % FALTAN LOS GIROS PARA ROTAR LA SECCION
    UG(2:6:end) = RotsAux(:,1) ;
    UG(4:6:end) = RotsAux(:,2) ;
    UG(6:6:end) = RotsAux(:,3) ;
  end
 


  if secc(1) == 25
  
    [ Conecvtk, Nodesvtk, vtkDispMat ] = vtkQuadHexa( NodesAux, ConecAux, secc, UG, Local2GlobalMatsDef, dispMat ) ;

  elseif secc(1) == 12     

    [ Conecvtk, Nodesvtk, vtkDispMat ] = vtkHexa    ( NodesAux, ConecAux, secc, UG, Local2GlobalMatsDef, dispMat ) ;

  else
    dispMat = reshape(UG(1:2:end),3,nnodes)' ;
    NodesDef = Nodes + [UG(1:6:end) UG(3:6:end) UG(5:6:end)] ;
    Nodesvtk = NodesDef ;
    Conecvtk = [ Conec(:,7) Conec(:,1:4) ] ;
    vtkDispMat = dispMat ;
    vtkNormalForceMat = [ vtkNormalForceMat ; normalForce ] ;
  end

end
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------


% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
currdir = pwd ;
lw  = 2   ; ms  = 5.5 ;
lw2 = 3.2 ; ms2 = 23 ;
plotfontsize = 22 ;

nelems = size(Conec) ;
nnodes = size(Nodes) ;
ndofpnode = 6 ;

% -----------------------------------------

for indplot = 1 : length( timesPlotsVec ) ;

  Utplot = matUts ( :, timesPlotsVec( indplot) ) ;

  [vtkConec, vtkNodesDef, vtkDispMat, vtkNormalForceMat ] = vtkConecNodes ( Nodes, Conec, sectPar , Utplot, nonLinearAnalysisBoolean, ...
                                                                            dynamicAnalysisBoolean, coordsElemsMat, matNts(:, timesPlotsVec(indplot)) ) ;

  filename = [ outputdir problemName '_' sprintf('%04i',indplot) '.vtk'] ;
	
  % Scalars vals
  svm             = [] ;
  vecSigI         = [] ;
	vecSigII        = [] ;
	vecSigIII       = [] ;
  normalForceMat  = [] ;
  
  % Vectors vals
  vecI    = [] ;
	vecII   = [] ;
	vecIII  = [] ;
  
  cellPointData   = cell ;
	cellCellData    = cell ;
  cellTensorData  = cell ;	
  
  
	if size(matNts,1) > 0
    cellCellData{1,1} = 'SCALARS' ; cellCellData{1,2} = 'Normal_Force' ; cellCellData{1,3} = vtkNormalForceMat ;
	end
    
  if size(cellStress,1) > 0 && size(cellStress{3},1) > 0
    stressMat = cellStress{3} ;
    sxx = stressMat(:,1);
    syy = stressMat(:,2);
    szz = stressMat(:,3);
    tyz = stressMat(:,4);
    txz = stressMat(:,5);
    txy = stressMat(:,6);
    svm = sqrt( 0.5*( (sxx-syy).^2 + (syy-szz).^2 + (szz-sxx).^2  + 6*(tyz.^2 + txz.^2 + txy.^2) ) ) ;
    for i = 1:nelems
      tensor = [ sxx(i) txy(i) txz(i) ; txy(i) syy(i) tyz(i) ; txz(i) tyz(i) szz(i) ] ;
      [vec,val] = eig(tensor) ;
      [sigI,indI] = max(diag(val)) ;
      [sigIII,indIII] = min(diag(val)) ;
      indII = setdiff([1 2 3], [indI indIII]) ;
      sigII = diag(val)(indII) ;
      % Sigma I
      vecSigI = [vecSigI ; sigI] ;
      vecI = [vecI ; vec(:,indI)'/(norm(vec(:,indI)'))*abs(sigI)] ;
      % Sigma II
      vecSigII = [ vecSigII ; sigII ] ;
      vecII = [ vecII ; vec(:,indII)'/(norm(vec(:,indII)'))*abs(sigII)] ;
      % Sigma III
      vecSigIII = [vecSigIII ; sigIII] ;
      vecIII = [vecIII ; vec(:,indIII)'/(norm(vec(:,indIII)'))*abs(sigIII)] ;
        
    end
    cellCellData{1,1} = 'SCALARS' ; cellCellData{1,2} = 'Von_Mises'   ; cellCellData{1,3} = svm ;
    cellCellData{2,1} = 'VECTORS' ; cellCellData{2,2} = 'vI'          ; cellCellData{2,3} = vecI ;
    cellCellData{3,1} = 'VECTORS' ; cellCellData{3,2} = 'vII'         ; cellCellData{3,3} = vecII ;
    cellCellData{4,1} = 'VECTORS' ; cellCellData{4,2} = 'vIII'        ; cellCellData{4,3} = vecII ;
    cellCellData{5,1} = 'SCALARS' ; cellCellData{5,2} = 'SigI'        ; cellCellData{5,3} = vecSigI ;
    cellCellData{6,1} = 'SCALARS' ; cellCellData{6,2} = 'SigII'       ; cellCellData{6,3} = vecSigII ;
    cellCellData{7,1} = 'SCALARS' ; cellCellData{7,2} = 'SigIII'      ; cellCellData{7,3} = vecSigIII ;
  end

  cellPointData = cell(1,3) ;
  cellPointData{1,1} = 'VECTORS' ; cellPointData{1,2} = 'Displacements' ; cellPointData{1,3} = vtkDispMat ;

  
  vtkWriter( filename, vtkNodesDef, vtkConec , cellPointData, cellCellData ) ;
  
end




