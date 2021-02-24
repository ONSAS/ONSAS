% ---------------------------
% Heat Transfer FEM prototype solver
% ---------------------------

function [Ts, NodesCoord, tiempos] = HeatFEM( ...
  timeIncr, Tfinal, ...
  materialParams, ...
  geometryParams, ...
  meshParams, ...
  hConv, diriDofs, robiDofs, Tamb, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initialTempFunc, internalHeatFunc, ...
  diriFacesAndVals, neumFacesAndVals, robiFacesAndVals  );

% ==============================================
% seteos previos

% variables auxiliares para futuro en ONSAS
close all

% material params
rho   = materialParams(1) ;  
cSpHe = materialParams(2) ;
kCond = materialParams(3) ;


%
geometryType = geometryParams(1) ;

if geometryType == 1 % unidimensional
  L = geometryParams(2); Area = geometryParams(3);
elseif geometryType == 2 % solid 
  Lx = geometryParams(2);
  Ly = geometryParams(3);
  Lz = geometryParams(4);
end

% mesh creation
if geometryType == 1
  ndivs = meshParams (1) ;
elseif geometryType == 2
  ndivs = meshParams (1:3) ;
  
  ndirifaces = size( diriFacesAndVals,1) ;
  nrobifaces = size( robiFacesAndVals,1) ;
  nneumfaces = size( neumFacesAndVals,1) ;
end

if geometryType == 1
  NodesCoord = linspace(0, L, ndivs(1)+1 )' ;
  Conec = [ (1:(ndivs(1)))' (2:(ndivs(1)+1))' ] ;
elseif geometryType == 2
  [NodesCoord, Conec] = createSolidRegularMesh( ndivs, Lx, Ly, Lz ) ;
end




% defino/renombro variables de entrada
dt = timeIncr ;

if exist('nt')==0, nTimes = Tfinal / dt + 1 ; end

tiempos = linspace( 0, Tfinal, nTimes )' ;

if isempty( nPlots ),
  nPlots = 4 ;
  plotBoolean = 1,

elseif nPlots == 0,
  plotBoolean = 0,

else
  plotBoolean = 1,
  if nPlots == inf,
    nPlots = nTimes ;
  end

end


nnodes = size(NodesCoord,1)
nelem  = size(Conec,1) ;



% ================================================
% dofs for 3D
if geometryType == 2

  nodesFaceOne = (1:(ndivs(1)+1):nnodes)' ; % x=0
  nodesFaceTwo = ((ndivs(1)+1):(ndivs(1)+1):nnodes)' ;  % x=Lx

  nodesFaceThree =  [ ] ; % y=0  
  for i=1:(ndivs(3)+1)
    nodesFaceThree = [ nodesFaceThree ; ...
     ( (i-1)*((ndivs(1)+1)*(ndivs(2)+1)) + (1:(ndivs(1)+1))' ) ] ;
  end  

  nodesFaceFour =  nodesFaceThree + ((ndivs(1))*(ndivs(2)+1)) ; % y=0  

  nodesFaceFive = (1:1:((ndivs(1)+1)*(ndivs(2)+1)) )' ; % z=0
  nodesFaceSix  = (ndivs(1)+1)*(ndivs(2)+1)*(ndivs(3)) + nodesFaceFive ; % z=Lx

  facesCell = { nodesFaceOne, nodesFaceTwo, ...
                nodesFaceThree, nodesFaceFour, ...
                nodesFaceFive, nodesFaceSix} ;
  
  triangAreaFactorsFaceOne   = patchTriangFactors( ndivs(2)+1, ndivs(3)+1, 1 ) * Lx/ndivs(2)*Ly/ndivs(3) * 0.5 ;
  triangAreaFactorsFaceTwo   = patchTriangFactors( ndivs(2)+1, ndivs(3)+1, 2 ) * Lx/ndivs(2)*Ly/ndivs(3) * 0.5 ;
  triangAreaFactorsFaceThree = patchTriangFactors( ndivs(1)+1, ndivs(3)+1, 3 ) * Lx/ndivs(1)*Ly/ndivs(3) * 0.5 ;
  triangAreaFactorsFaceFour  = patchTriangFactors( ndivs(1)+1, ndivs(3)+1, 4 ) * Lx/ndivs(1)*Ly/ndivs(3) * 0.5 ;
  triangAreaFactorsFaceFive  = patchTriangFactors( ndivs(1)+1, ndivs(2)+1, 5 ) * Lx/ndivs(1)*Ly/ndivs(2) * 0.5 ;
  triangAreaFactorsFaceSix   = patchTriangFactors( ndivs(1)+1, ndivs(2)+1, 6 ) * Lx/ndivs(1)*Ly/ndivs(2) * 0.5 ;
  
  cellTriangFactors = { triangAreaFactorsFaceOne , triangAreaFactorsFaceTwo , triangAreaFactorsFaceThree , ...
                        triangAreaFactorsFaceFour, triangAreaFactorsFaceFive, triangAreaFactorsFaceSix   };
  
  if ndirifaces > 0
    diriDofs = unique( addFacesDofs( diriFacesAndVals(:,1), facesCell ) ) ;
  end
  if nrobifaces >0
    robiDofs = unique( addFacesDofs( robiFacesAndVals(:,1), facesCell ) ) ;
  end

end

% neuman dofs
neumdofs = 1:nnodes ;
neumdofs( diriDofs ) = [] ;
% ==============================================


% ==============================================
% construccion matrices FEM

if geometryType == 1
  % calculo largo de cada elemento
  lelem  = L/nelem ;
  
  % local elemental diffussion equation
  Kdiffe = kCond * Area / lelem * [ 1 -1 ; -1 1 ] ;
  
  %~ MintEe = rho * cSpHe * Area * lelem / 6 * [ 2 1 ; 1 2 ] ;
  MintEe = rho * cSpHe * Area * lelem / 2 * [ 1 0 ; 0 1 ] ;
  
  % unitary volumetric external heat source
  bQhe = 0.5 * Area * lelem * [ 1 ; 1 ] ;

  % initial temperature
  T0 = feval( initialTempFunc, NodesCoord ) ;
elseif geometryType == 2

  % initial temperature
  T0 = feval( initialTempFunc, NodesCoord(:,1) ) ;
  
end


% ------------------------
% Ts matriz para almacenar temperatuas en cada tiempo
Ts = T0 ;

% matrices/vectors assembly
KdiffG = sparse( nnodes, nnodes ) ;
MintEG = sparse( nnodes, nnodes ) ;
MrobiG = sparse( nnodes, nnodes ) ;

QhG    = zeros( nnodes, 1      ) ;

for i = 1 : nelem

  nodeselem = Conec(i,:) ;
  elemDofs  = nodes2dofs ( nodeselem, 1 ) ;

  if geometryType == 2
    [ Kdiffe, MintEe, bQhe ] = elementHeatTetraSolid( NodesCoord( nodeselem,:) , materialParams ) ;
  end  

  KdiffG( elemDofs , elemDofs ) = ...
  KdiffG( elemDofs , elemDofs ) + Kdiffe  ;

  MintEG( elemDofs , elemDofs ) = ...
  MintEG( elemDofs , elemDofs ) + MintEe  ;

  QhG   ( elemDofs            ) = ...
  QhG   ( elemDofs            ) + bQhe  ;

end

if  ~isempty( robiDofs )
  if geometryType == 1
    MrobiG ( robiDofs, robiDofs ) = hConv ;
  elseif geometryType == 2
    if nrobifaces > 0
      for k=1:nrobifaces
        currFace = robiFacesAndVals (k,1) ;
        MrobiG ( facesCell{currFace}, facesCell{currFace} ) = MrobiG ( facesCell{currFace}, facesCell{currFace} ) ...
                                                            + diag( cellTriangFactors{currFace} ) * robiFacesAndVals(k,2) ;
      end
    end
  end
  KdiffG = KdiffG + MrobiG ;
end


% ------------------------

CDD = MintEG( diriDofs, diriDofs ) ;
CND = MintEG( neumdofs, diriDofs ) ;
CNN = MintEG( neumdofs, neumdofs ) ;

KdiffGNN = KdiffG(neumdofs, neumdofs ) ;
KdiffGND = KdiffG(neumdofs, diriDofs ) ;

Matrix = ( KdiffGNN * dt + CNN )  ;

qext = zeros( nnodes, 1 ) ;

% input fluxes
if ~isempty( qInpLeft  ), qext(   1) = qInpLeft ; end
if ~isempty( qInpRight ), qext( end) = qInpRight; end

if ~isempty( diriDofs )
  qext( neumdofs ) = qext( neumdofs ) - KdiffGND * ones( length( diriDofs ), 1 ) * Tdiri ;
end

if ~isempty( robiDofs )
  if geometryType == 1
    qext( robiDofs ) = qext( robiDofs ) + hConv * Tamb ;
  elseif geometryType == 2
  
    if nrobifaces > 0
      for k=1:nrobifaces
        currFace = robiFacesAndVals(k,1) ;
        qext ( facesCell{currFace} ) = qext ( facesCell{currFace} ) + cellTriangFactors{currFace} * robiFacesAndVals(k,2) * Tamb ;
      end
    end
  end
end

qext = qext ;
% ==============================================


% =========================================
% simulation

if plotBoolean
  indsPlotStep = round( nTimes / nPlots ) ;
  if length( ndivs)==1
    figure, hold on, grid on
  end
end

fext = zeros( size(Ts) ) ;

% ------------------------
% time loop

for ind = 1:nTimes %ind es el indice de tiempo que se esta hallando
  t = dt*(ind-1) ;
  fprintf('ind: %4i  time: %15.4e\n',ind, t);
  
  if ~isempty( internalHeatFunc )
    fext ( neumdofs, ind ) = qext ( neumdofs    ) + QhG(neumdofs) * feval( internalHeatFunc, t ) ;
  else
    fext ( neumdofs, ind ) = qext ( neumdofs    ) ;
  end
  

  if ind > 1
    f = (    fext( neumdofs, ind ) * dt   ...
         + MintEG( neumdofs, : ) * Ts( :, ind-1 )
        );

    Tip1 = Matrix \ f ;
    Ts( neumdofs, ind ) = Tip1 ;

    if ~isempty( diriDofs ),
      Ts( diriDofs, ind ) = Tdiri   ;
    end
    
  end % if i not zero

  % --- plots ---
  if plotBoolean
    if ind == 1 || mod( ind, indsPlotStep ) == 0
      plotTemp( NodesCoord, Conec, ind, Ts(:,ind), geometryType )
    end
  end
  % ---------------

end % for times
% ------------------------

if plotBoolean && length(ndivs)==1
  print( [ 'pngs/' problemName '.png'],'-dpng')
end

M = sparse( nnodes, nnodes ) ;
K = KdiffGNN ;
C = CNN ;
xs = NodesCoord(:,1) ;


us    = Ts ;
udots = [] ;

save('-mat', ['mats/' problemName '.mat'], 'K', 'C', 'M', 'fext', 'timeIncr', 'us', 'udots', 'xs')





function factors = patchTriangFactors( nnod1, nnod2, face );

if face == 1 || face == 2 || face == 5 || face == 6
  booleanSense12 = true;
else
  booleanSense12 = false;
end

factors = zeros( nnod1* nnod2 ,1 ) ;

for j=1:nnod2
  inds = (j-1)*(nnod1) + (1:(nnod1))' ;
  
  if j == 1
    factors( inds(2:(end-1)) ) = 1; 
    if booleanSense12
      factors( inds(1)   ) = 1/3 ; 
      factors( inds(end) ) = 2/3 ; 
    end

  elseif j == nnod2
    factors( inds(2:(end-1)) ) = 1; 
    if booleanSense12
      factors( inds(1)   ) = 2/3 ; 
      factors( inds(end) ) = 1/3 ; 
    end
    
  else
    factors( inds(2:(end-1)) ) = 2; 
    factors( inds(1)   ) = 1 ; 
    factors( inds(end) ) = 1 ;   
  end
  
end



function [dofsVec] = addFacesDofs( faces, facesCell )
  
  dofsVec = [];
  if length(faces) > 0
    for i=1:(length(faces))
      dofsVec = [ dofsVec; facesCell{faces(i)} ] ;
    end
  end

