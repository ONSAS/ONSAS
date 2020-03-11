% ------------------------------------------------------------------------------
function [Conecvtk, Nodesvtk, vtkDispMat, vtkNormalForceMat ] = vtkConecNodes ( Nodes, Conec, indexesElems, secc, UG, nonLinearAnalysisBoolean, dynamicAnalysisBoolean, coordsElem, normalForce )
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
        titax = [] ; titay = [] ; titaz = [] ;
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
    % Replaces nodes x, y, z rots
    if size(RotsAux) > 0
      UG(2:6:end) = RotsAux(:,1) ;
      UG(4:6:end) = RotsAux(:,2) ;
      UG(6:6:end) = RotsAux(:,3) ;
    end
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
    for i = 1:nelems
      m = indexesElems(i) ;
      if Conec(i,7) == 1 || Conec(i,7) == 2
        vtkNormalForceMat = [ vtkNormalForceMat ; normalForce(m) ] ;
      elseif Conec(i,7) == 3 || Conec(i,7) == 4
        vtkNormalForceMat = [ vtkNormalForceMat ; 0 ] ;
      end  
    end    
  end


end % end function declaration
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
