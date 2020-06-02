% ------------------------------------------------------------------------------
% Test problem: Plate submitted to distributed load on surface
% ------------------------------------------------------------------------------

inputONSASversion = '0.1.9';

problemName = 'Plate' ;

ndofpnode = 6 ;

% Material properties
E   = 3004160 ; 
nu  = 0.3 ;

hyperElasParams = cell(1,1) ;  
hyperElasParams{1} = [ 1 E nu ] ;

% Sections
L   = [ 4 4 ]' ;  % L = [ length in dimension "x", length dimension "y"]
t   = 0.1   ;     % thickness
secGeomProps = [ 1 1 1 1 ] ;

% Mesh
nx  = 34 ;                    % # divitions in direction "x"

ny      = ceil(L(2)/L(1)*nx) ;  % # divitions in direction "y"
nel     = [ nx ny ]' ;
nnos    = nel + 1 ;  % # of nodes in each direction
nnostot = nnos(1)*nnos(2) ;  % # of total nodes

% Nodes
lins1   = linspace( 0 , L(1) , nnos(1) )';
lins2   = linspace( 0 , L(2) , nnos(2) )';

Nodes = [ ] ;  
for i = 1:nnos(2) 
  Nodes( ( (nnos(1) * (i-1) + 1) : nnos(1)*i ) ,:) = [ lins1 lins2(i)*ones(nnos(1),1) zeros(nnos(1),1) ] ;
end

% Conectivity matrix 
Conec = [ ] ; 
for j = 1:nel(2)
  for i = 1:nel(1)
    intri = (i-1)+1+(j-1)*nel(1) ;
    Conec( intri , : ) = [(j-1)*nnos(1)+i    (j-1)*nnos(1)+i+1   j*nnos(1)+i+1   j*nnos(1)+i] ;
  end
end

Conec = [ Conec ones(size(Conec,1),1) zeros(size(Conec,1),1) ones(size(Conec,1),1)*4 ] ;

% Loads 

q = -1 ;

nnodes = size(Nodes,1) ;
nelems = size(Conec,1) ;
midnode = floor(nnodes/2) +1 ;
a = L(1) / (nel(1) * 2) ;
b = L(2) / (nel(2) * 2) ;
%~ qelemaux = 4*q*a*b* [1/4   a/12   b/12   1/4   -a/12    b/12   1/4   -a/12   -b/12   1/4   a/12   -b/12 ] ;
%~ qelem = 4*q*a*b* [ 0 a/12 0 b/12 1/4 0 0 -a/12 0 b/12 1/4 0 0 -a/12 0 -b/12 1/4 0 0 a/12 0 -b/12 1/4 0 ]' ;
%~ qtot  = zeros(ndofpnode*nnodes,1) ;

%~ for i = 1:nelems
  %~ if Conec(i,7) == 4
  %~ nodeselem = Conec(i,1:4) ;
  %~ elemdofs = nodes2dofs(nodeselem,ndofpnode) ;
  %~ qtot (elemdofs) = qtot(elemdofs) + qelem ;
  %~ end
%~ end

%~ nodalConstantLoads = zeros(nnodes,7) ;
%~ for i = 1:nnodes
%~ dofs = nodes2dofs(i,ndofpnode) ;
%~ nodalConstantLoads(i,:) = [i qtot(dofs)' ] ;
%~ end

% Boundary Conditions

NodesSouth  = ( 1:nnos(1) )' ;
NodesEast   = ( nnos(1):nnos(1):nnodes )' ;
NodesNorth  = ( nnos(1)*ny+1:nnostot    )' ;
NodesWest   = ( 1:nnos(1):nnos(1)*ny+1  )' ;

% fixed degrees of freedom 

% Boundary conditions:  0 free 1 displacement fixed 2 displacement and rotation fixed
CondSouth = 1 ;
CondNorth = 1 ;
CondWest  = 1 ;
CondEast  = 1 ;

nodalSprings = [] ;


for i = 1:nnodes

  if ismember(i,NodesSouth)
    if CondSouth == 0 
  
    elseif CondSouth == 1
      nodalSprings(i,:) = [i inf 0 inf 0 inf 0 ] ;
    else 
      nodalSprings(i,:) = [i ones(1,6)*inf] ;
    end
    
  elseif ismember(i,NodesEast)
    if CondEast == 0
    
    elseif CondEast == 1
      nodalSprings(i,:) = [i inf 0 inf 0 inf 0 ] ;
    else
      nodalSprings(i,:) = [i ones(1,6)*inf] ;
    end
  elseif ismember(i,NodesWest)
    if CondWest == 0
    
    elseif CondWest == 1
      nodalSprings(i,:) = [i inf 0 inf 0 inf 0 ] ;
    else 
      nodalSprings(i,:) = [i ones(1,6)*inf] ;
    end
  elseif ismember(i,NodesNorth) 
    if CondNorth == 0
  
    elseif CondNorth == 1 
      nodalSprings(i,:) = [i inf 0 inf 0 inf 0 ] ;
    else 
      nodalSprings(i,:) = [i ones(1,6)*inf] ;
    end  
  end

end

null = find(nodalSprings(:,1)==0) ;
nodalSprings(null,:) = [] ;

% Analysis parameters

nonLinearAnalysisBoolean = 0 ; 
printflag = 2 ;
linearDeformedScaleFactor = 10.0 ;

% Plot options

plotParamsVector = [ 3 ] ;
printflag = 2 ;

% Analytic sol flag

analyticSolFlag = 3 ;

alpha1 = pi/2 ; 
alpha3 = 3*pi/2 ;

m1 = 1 ;
fu1 = ((-1)^((m1-1)/2) * (alpha1 * tanh(alpha1) + 2) / (m1^5 * 2*cosh(alpha1)) ) ;
m3 = 3 ;
fu3 = ((-1)^((m3-1)/2) * (alpha3 * tanh(alpha3) + 2) / (m3^5 * 2*cosh(alpha3)) ) ;

D = E * t^3 / ( 12 * (1-nu^2) ) ;

w1 = 5/384 ;
const =(4/pi^5) ;
w2 = -const*(fu1+fu3) ;

analytSol = (w1 + w2)*q*L(1)^4/D ;
analyticSolDofs = [ midnode*ndofpnode-1 ] ;
analyticCheckTolerance = 1e-3 ;

nelems = size(Conec,1) ;
nelemsVec = [1:nelems]' ;
unifLoad = [ nelemsVec ones(nelems,1) zeros(nelems,2) ones(nelems,1)*q ] ;



