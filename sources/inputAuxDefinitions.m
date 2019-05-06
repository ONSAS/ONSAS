
% ---------------------------------------------------

ndofpnode = 6;


%~ strucsize = max( [ max(Nodes(:,1))-min(Nodes(:,1)) , ...
                   %~ max(Nodes(:,2))-min(Nodes(:,2)) , ...
                   %~ max(Nodes(:,3))-min(Nodes(:,3)) ] ) ; 

strucsize = strucSize(Nodes) ;

absoluteimperfection = imperfactor * strucsize ;

if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
  controlDofInfo = [ ] ;
else 
  % control dof info
  controlDof       = nodes2dofs( controlDofInfo(1), 6)(controlDofInfo(2)) ;
  controlDofFactor = controlDofInfo(3) ; 
end


releasesDofs = [];
for i=1:nelems
  if Conec(i,7)==1
    releasesDofs = [ releasesDofs; nodes2dofs( Conec(i,1:2),ndofpnode)(2:2:end) ];
    Releases = [ Releases; i ones(1,4) ] ;
  end
end

releasesDofs = unique( releasesDofs);

coordsElemsMat = zeros(nelems,2*ndofpnode) ;
for i=1:nelems
  % obtains nodes and dofs of element
  nodeselem = Conec(i,1:2)' ;
  dofselem  = nodes2dofs( nodeselem , ndofpnode ) ;
  coordsElemsMat( i, (1:2:11) ) = [ Nodes( nodeselem(1),:)  Nodes( nodeselem(2),:) ] ;
end


hyperElasParamsMat = [] ;
for i=1:size(hyperElasParams,1)
  hyperElasParamsMat (i,1:length(hyperElasParams{i})) = hyperElasParams{i} ;
end



% ---------------- load vectors assembly -----------------------
variableFext = zeros( ndofpnode*nnodes , 1 );
constantFext = zeros( ndofpnode*nnodes , 1 );

if exist( 'nodalVariableLoads' ) ~= 0
  for i=1:size(nodalVariableLoads,1)
    aux = nodes2dofs ( nodalVariableLoads(i,1), ndofpnode ) ;
    variableFext( aux ) = variableFext( aux ) + nodalVariableLoads(i,2:7)' ;
  end
end

if exist( 'nodalConstantLoads' ) ~= 0
  for i=1:size(nodalConstantLoads,1)
    aux = nodes2dofs ( nodalConstantLoads(i,1), ndofpnode ) ;
    constantFext( aux ) = constantFext( aux ) + nodalConstantLoads(i,2:7)' ;
  end
end


% ------------------------------------------------------------
if exist( 'nodalConstantLoads' ) ~= 0 || exist( 'nodalVariableLoads' ) ~= 0
  [maxNorm2F, visualloadfactor] = visualLoadFac(strucsize, variableFext, constantFext, nnodes) ;
else
  visualloadfactor = 1 ;
end

% -------------------- indexes computation --------------------
indexesElems = zeros(nelems,1) ;
ntruss = 0 ;
nbeam = 0 ;
ntet = 0 ;
nplate = 0 ;
trussElem = [] ;
beamElem = [] ;
beamNodes = [] ;
tetElem = [] ;
plateElem = [] ;

for i = 1:nelems
  if Conec(i,7) == 1
    nodeselem = Conec(i,1:2) ;
    ntruss = ntruss + 1 ;
    indexesElems(i) = ntruss ;
    trussElem = [ trussElem ; i ] ;
  elseif Conec(i,7) == 2
    nodeselem = Conec(i,1:2) ;
    nbeam = nbeam + 1 ;
    indexesElems(i) = nbeam ;
    beamElem = [ beamElem ; i ] ;
    beamNodes = [ beamNodes ; nodeselem' ] ;
  elseif Conec(i,7) == 3
    ntet = ntet + 1 ;
    indexesElems(i) = ntet ;
    tetElem = [ tetElem ; i ] ;
  elseif Conec(i,7) == 4
    nplate = nplate + 1 ;
    indexesElems(i) = nplate ;
    plateElem  = [ plateElem ; i ] ;
  end
end

% Output parameters

cellStress = [] ;
matNts = [] ;
matUts = [] ;
% ------------------------------------------------------------

