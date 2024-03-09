function [ Waket1p, WakeQSt1p, WakeINTt1p ] = BEMcomputeInducedVel( BCsData, modelProperties, modelCurrSol, nextTime )

neumdofs = BCsData.neumDofs ;
Conec    = modelProperties.Conec;  nElems   = size(Conec, 1) ;
elements = modelProperties.elements; 
Nodes    = modelProperties.Nodes;  nNodes   = size(Nodes, 1) ;

timeIndexWake  = modelCurrSol.timeIndex + 1 ;
currTime       = nextTime ;

DWMbool    = modelProperties.analysisSettings.dynFlowModel;
WakeQSt1p  = zeros( 3*nNodes, 1 ) ;
WakeINTt1p = zeros( 3*nNodes, 1 ) ;
Waket1p    = zeros( 3*nNodes, 1 ) ;

for elem = 1:nElems
    mebVec     = Conec( elem, 1:3) ;
    % extract element properties
    elemType   = elements( mebVec( 2 ) ).elemType     ;
    modelName  = modelProperties.materials.modelName  ;

    % obtain element info
    [numNodes, nodalDofsEntries] = elementTypeDofs( elemType ) ;

    % obtains nodes and dofs of element
    nodeselem   = Conec( elem, (3+1):(3+numNodes) )' ;
    dofselem    = nodes2dofs( nodeselem , 6 )   ;
    dofsWake    = nodes2dofs( nodeselem , 3 )   ;

    % construct vector of degrees of freedom of element
    auxA = repmat( nodalDofsEntries, length(dofselem)/6,1 )  ;
    auxB = repelem( (0:6:length(dofselem)-1)',length(nodalDofsEntries),1) ;
    dofselemRed = dofselem( auxA+auxB )   ;

    %md elemDisps contains the displacements corresponding to the dofs of the element
    elemDisps       = Utp1( dofselemRed )   ;
    dotdispsElem    = Udottp1(dofselemRed ) ;

    elemNodesxyzRefCoords  = reshape( Nodes( nodeselem, : )', 1, 3*numNodes ) ;
    
    if strcmp( elemType, 'frame')
        if strcmp( modelName, 'elastic-rotEngStr') || strcmp(modelName, 'elastic-linear')
            if ~isempty(modelProperties.elements.BEMparams) && BEMbool
                elemWaket       = Wake( dofsWake    ) ;
                elemWakeQSt     = WakeQS( dofsWake  ) ;
                elemWakeINTt    = WakeINT( dofsWake ) ;

                [ Waket1p(dofselemRed), WakeQSt1p(dofselemRed), WakeINTt1p(dofselemRed) ] = uBEMupdateInducedVelocity(modelProperties.analysisSettings, elemNodesxyzRefCoords   , ...
                                                                                                    elemWaket, elemWakeINTt, elemWakeQSt                   , ...
                                                                                                    elemDisps, dotdispsElem,                 ...
                                                                                                    modelProperties.elements.BEMparams ,     ...
                                                                                                    modelProperties.elements.airFoilPolars,  ...
                                                                                                    modelProperties.elements.dynStallParams,  ...
                                                                                                    currTime, DWMbool);

            end
        end
    end
end