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


% ==============================================================================

%Script for Linear Analysis.

if plotParamsVector(1)>0
  fprintf('  - performing linear analysis ... ');
end

% ==============================================================================
% ---------------------    Geometry reading    ------------------------
tic ;
ndofpnode = 6 ;

ElemLengths = zeros(nbeam+ntruss,1) ;
tetVol      = zeros(ntet,1) ;

Local2GlobalMats = cell(nbeam+ntruss,1) ;

eyetres      = eye(3)     ;
eyevoig      = zeros(6,1) ;
eyevoig(1:3) = 1.0        ;

BMat    = cell(ntet,1) ;
ConsMat = cell(ntet,1) ;

Lx = zeros(nplate,1) ;
Ly = zeros(nplate,1) ;

for i = 1:nelems
  m =  indexesElems(i) ;    
  
  if Conec(i,7) == 1 || Conec(i,7) == 2                                              
		
    [ ElemLengths(m) Local2GlobalMats{m} ] = beamParameters(Nodes(Conec(i,1:2),:)) ; 
  
  elseif Conec(i,7) == 3
    
    nodeselem = Conec(i,1:4) ;
    dofselem  =  nodes2dofs ( nodeselem , ndofpnode   ) ;
    
    dofstet = dofselem(1:2:length(dofselem)) ;
    
    tetcoordmat        = zeros(3,4) ;
    tetcoordmat(1,1:4) = Nodes( nodeselem , 1 ) ;
    tetcoordmat(2,1:4) = Nodes( nodeselem , 2 ) ;
    tetcoordmat(3,1:4) = Nodes( nodeselem , 3 ) ;
    
    [ deriv , vol ] =  DerivFun( tetcoordmat ) ;
    
    if vol<0, i, error('Element with negative volume, check connectivity.'), end
    
    tetVol(m) = vol ;
    
    BMat{m} = BMats ( deriv ) ;
    
    E  = hyperElasParams{ Conec(i,5)}(1+1) ;  
    nu = hyperElasParams{ Conec(i,5)}(1+2) ; 
    
    mu    = E / (2.0 * ( 1.0 + nu ) ) ;
    Bulk  = E / (3.0 * ( 1.0 - 2.0 * nu ) ) ;
    
    ConsMataux = zeros(6,6) ;
    ConsMataux (1:3,1:3) = Bulk  + 2* mu * ( eyetres - 1.0/3.0 ) ;
    ConsMataux (4:6,4:6) =            mu *   eyetres ;
    
    ConsMat{m} = ConsMataux ;
  
  elseif Conec(i,7) == 4 
    
    nodeselem = Conec(i,1:4) ;
    nod1  = nodeselem(1) ;
    nod2  = nodeselem(2) ;
    nod4  = nodeselem(4) ;
    Lx(m) = abs(Nodes(nod1,1) - Nodes(nod2,1)) ;
    Ly(m) = abs(Nodes(nod1,2) - Nodes(nod4,2)) ;
  
  end  
end
tGeomReading = toc ;
% ==============================================================================


% ==============================================================================
% ---------------------    Stiffness Matrix Assembly    ------------------------

if plotParamsVector(1)>0
  fprintf('Assembling stiffness matrix...\n');
end
tic ;
KG      = sparse( ndofpnode*nnodes, ndofpnode*nnodes ) ;
ElemKGs = cell(nbeam+ntruss,1) ;                                        

elemReleases = zeros(nbeam+ntruss,4) ;
gdlReleases = [] ;  
                    
if exist('Releases') ~= 0
  if size(Releases,1)>0
    for i=1:size(Releases,1)
      elemReleases( Releases(i,1), : ) = Releases(i, 2:5 ) ; 
    end
  end
end


for i = 1:nelems
  m = indexesElems(i) ;
  if mod(m,round(nelems/20)) == 0 ,
    if plotParamsVector(1)>0, fprintf('=') ; end
  end

  if Conec(i,7) == 1 || Conec(i,7) == 2                                 
 
    % material properties
    E  = hyperElasParams{ Conec(i,5)}(1+1) ;  
    nu = hyperElasParams{ Conec(i,5)}(1+2) ;  
    % ------------------------
  
    % section properties
    A  = secGeomProps ( Conec(i,6), 1 ) ;  
    Iy = secGeomProps ( Conec(i,6), 2 ) ;  
    Iz = secGeomProps ( Conec(i,6), 3 ) ;  
    J  = secGeomProps ( Conec(i,6), 4 ) ;  
    l  = ElemLengths(m) ;
    
    % -------------------------
    % sets the nodes of the element and the corresponding dofs
    nodi = Conec(i,1) ;    nodj = Conec(i,2) ;  
    elemdofs = nodes2dofs ( [ nodi nodj ]' , ndofpnode ) ;
    % --------------------------------
		
    R = RotationMatrix ( ndofpnode, Local2GlobalMats{m} ) ;
		
    KGelem = linearStiffMatBeam3D(E, nu, A, Iy, Iz, J, l, elemReleases(m,:), R) ;
    ElemKGs{i} = KGelem;

  elseif Conec(i,7) == 3
     
    nodeselem = Conec(i,1:4) ;
    elemdofs  =  nodes2dofs ( nodeselem , ndofpnode ) ;
   
    Kml    = BMat{m}' * ConsMat{m} * BMat{m} * tetVol(m) ;
    KGelem = zeros(24,24) ;  
    KGelem([1:2:end], [1:2:end])  =  Kml ;

  elseif Conec(i,7) == 4 % plates
  
    nodeselem = Conec(i,1:4) ;
    elemdofs = nodes2dofs( nodeselem, ndofpnode ) ;
    lx = Lx(m) ;
    ly = Ly(m) ;
    
    Kelem = linearStiffMatPlate3D(E, nu, t, lx, ly) ;  
    KGelem = zeros(24,24) ;
    dofs = [5,2,4,5+ndofpnode,2+ndofpnode,4+ndofpnode,5+2*ndofpnode,2+2*ndofpnode,4+2*ndofpnode,5+3*ndofpnode,2+3*ndofpnode,4+3*ndofpnode] ;
    KGelem(dofs,dofs) = Kelem ;
   
  end

  KG( elemdofs     , elemdofs     ) = KG( elemdofs , elemdofs )  +  KGelem ;

end % endfor elems


% Fixed displacements and Springs 

KS = sparse( ndofpnode*nnodes, ndofpnode*nnodes ) ;

fixeddofs        = [ ] ;
fixeddofsR       = [ ] ;
fixeddofsSprings = [ ] ;


for i = 1:size(nodalSprings,1)
  aux = nodes2dofs ( nodalSprings (i,1) , ndofpnode ) ;

  for k = 1:ndofpnode
    %
    if nodalSprings(i,k+1) == inf 
			fixeddofsR = [ fixeddofsR; aux(k) ] ;
    elseif nodalSprings(i,k+1) ~= 0
      fixeddofsSprings = [ fixeddofsSprings ; aux(k)] ;
      KS(aux(k), aux(k) ) = nodalSprings(i,k+1) ;
    end
    
  end
end

% add springs stiffness matrix
KG = KG + KS ;
tStiffMatrix = toc ;
% ==============================================================================




% ==============================================================================
% ------------------------    Load Vector Assembly    --------------------------
tic ;
% Parameters 
if length(userLoadsFilename) > 0 || norm(variableFext) > 0
	deltaT = numericalMethodParams(2) ;
	tf = numericalMethodParams(3) ;
	indexTimeSol = numericalMethodParams(4) ;
	nTimeSteps = tf/deltaT + 1 ;
else
	numericalMethodParams = [] ;
	deltaT = 1 ;
	tf = 1 ;
	indexTimeSol = 1 ;
end

nodalConstUnif = zeros(ndofpnode*nnodes,1) ;
selfWeight = zeros(ndofpnode*nnodes,1) ;
nodalPlate = zeros(ndofpnode*nnodes,1) ;
elemUnifAux = zeros(nelems,ndofpnode*2) ;
elemSelfWAux = zeros(nelems,ndofpnode*2) ;
elemUnifLoc = zeros(nelems, ndofpnode*2) ;
elemSelfWLoc = zeros(nelems, ndofpnode*2) ;

% Nodal constant loads
[ nodalConstUnif, selfWeight, nodalPlate, ...
  elemUnifLoc, elemSelfWLoc, elemUnifAux, elemSelfWAux] = assemblyUniform(unifLoad, selfWeightBoolean, hyperElasParams, Conec, Nodes, indexesElems, ElemLengths, ...
                                                                          Local2GlobalMats, elemReleases, secGeomProps, Lx, Ly)  ;


% Assembly
constantFext = constantFext + nodalConstUnif + nodalPlate ; % Computes constant fext as the sum of the constant fext plus distributed loads

if length(userLoadsFilename) > 0
  matUG = zeros( ndofpnode*nnodes, nTimeSteps ) ; 
  Reactions = zeros( nnodes*ndofpnode, nTimeSteps ) ;
  matFs = [] ;
  for currTime = 1:nTimeSteps
    currLoadFactor = loadFactorsFunc(currTime) ;
    matFs = [ matFs (selfWeight + constantFext + currLoadFactor*variableFext + feval(userLoadsFilename , currTime-1 )) ] ;
  end

else 
% Variable forces
  if norm(variableFext) > 0
    if norm(selfWeight) > 0 && norm(constantFext) > 0 
      matFs = [ selfWeight constantFext variableFext ] ;
      matUG = zeros( ndofpnode*nnodes, 3 ) ;
      Reactions = zeros(nnodes*ndofpnode, 3) ;
    elseif norm(selfWeight) == 0 && norm(constantFext) > 0
      matFs = [ constantFext variableFext ] ;
      matUG = zeros( ndofpnode*nnodes, 2 ) ;
      Reactions = zeros(nnodes*ndofpnode, 2) ;
    elseif norm(selfWeight) > 0 && norm(constantFext) == 0
      matFs = [ selfWeight variableFext ] ;
      matUG = zeros( ndofpnode*nnodes, 2 ) ;
      Reactions = zeros(nnodes*ndofpnode, 2) ;
    else
      matFs = [ variableFext ] ;
      matUG = zeros( ndofpnode*nnodes, 1 ) ;
      Reactions = zeros(nnodes*ndofpnode, 1) ;  
    end  
    loadFactors = [] ;
    for currTime = 1:nTimeSteps
      loadFactors = [loadFactors loadFactorsFunc(currTime)] ;
    end
% Constant forces    
  else
    if norm(selfWeight) > 0 && norm(constantFext) > 0
      matFs = [ selfWeight constantFext ] ;
      matUG = zeros( ndofpnode*nnodes, 2 ) ;
      Reactions = zeros(nnodes*ndofpnode, 2) ;
      nTimeSteps = 2 ;
    elseif norm(selfWeight) == 0 && norm(constantFext) > 0
      matFs = [ constantFext ] ;
      matUG = zeros( ndofpnode*nnodes, 1 ) ;
      Reactions = zeros(nnodes*ndofpnode, 1) ;
      nTimeSteps = 1 ;
    else
      matFs = [ selfWeight ] ;
      matUG = zeros( ndofpnode*nnodes, 1 ) ;
      Reactions = zeros(nnodes*ndofpnode, 1) ;  
			nTimeSteps = 1 ;
    end
  end

end %end first if
tLoadsAssembly = toc ;
% ==============================================================================
% --------------------------    System Resolution   ----------------------------

% Non-zero prescribed displacements 
tic
fixeddofsD = [ ] ;
for i = 1:size(prescribedDisps,1)
  aux = nodes2dofs ( prescribedDisps (i,1), ndofpnode ) ;
  fixeddofsD = [ fixeddofsD ; aux(prescribedDisps(i,2) , 1 ) ] ;
  matUG ( aux(prescribedDisps(i,2) , : ) ) = prescribedDisps ( i, 3 ) ;   
end


% Dofs array
trussdofs = [];
tetdofs = [] ;
platedofs = [] ;

for i = 1:nelems
  if Conec(i,7) == 1
		nodeselem = Conec(i,1:2) ;
		elemdofs = nodes2dofs(nodeselem,ndofpnode) ;
		if dim == 1
			trussdofs = unique([ trussdofs ; elemdofs([2 3 4 5 6]) ; elemdofs([ 8 9 10 11 12]) ]) ;
		elseif dim == 2
			trussdofs = unique([ trussdofs ; elemdofs([2 3 4 6]) ; elemdofs([ 8 9 10 12]) ]) ;
		elseif dim == 3
			if ~ismember(nodeselem(1), beamNodes) && ~ismember(nodeselem(2), beamNodes)				
				trussdofs = unique([ trussdofs ; elemdofs([2 4 6]) ; elemdofs([ 8 10 12]) ]) ;
			elseif ~ismember(nodeselem(1), beamNodes) 
				trussdofs = unique([ trussdofs ; elemdofs([2 4 6]) ]) ;
			elseif ~ismember(nodeselem(2), beamNodes)
				trussdofs = unique([ trussdofs ; elemdofs([8 10 12]) ]) ;
			end
		end
		
  elseif Conec(i,7) == 3
    nodeselem = Conec(i,1:4) ;
    elemdofs = nodes2dofs(nodeselem,ndofpnode) ;
    tetdofs = [ tetdofs ; elemdofs(2:2:end) ] ;
  elseif Conec(i,7) == 4
    nodeselem = Conec(i,1:4) ;
    for j = 1:4
      if ismember(nodeselem(j), beamNodes)
        platedofs = platedofs ;
      else
        platedofs = [ platedofs ; nodeselem(j)*ndofpnode-5 ; nodeselem(j)*ndofpnode-3 ; nodeselem(j)*ndofpnode ] ;
      end
    end
  end
  
end


fixeddofs = unique([fixeddofsR ; fixeddofsD ; trussdofs ; tetdofs ; platedofs ]) ;
notfixeddofs             = 1:( ndofpnode*nnodes ) ;
notfixeddofs ( fixeddofs ) = [ ] ;
% Stiffness dofs

Kliblib = KG ( notfixeddofs , notfixeddofs ) ;
Klibcon = KG ( notfixeddofs , fixeddofs ) ; 
Kconlib = KG ( fixeddofs , notfixeddofs ) ;
Kconcon = KG ( fixeddofs , fixeddofs ) ; 



Flib = matFs ( notfixeddofs, : ) ; 

if length(fixeddofsSprings) == 0 % Check if there is any spring support (not rigid)
  Ucon = matUG ( fixeddofs, : ) ; % If there is any spring support with k<inf Ucon consider that node
else
  Ucon = [ ] ; % If it is not, Ucon is computed as an empty array
end

if length(Ucon) == 0
  Ulib = Kliblib \ Flib ;
  Fcon = matFs (fixeddofs, :) ;
else
  Ulib = Kliblib \ ( Flib - Klibcon*Ucon ) ; 
  Fcon = Kconlib*Ulib + Kconcon*Ucon ; 
end



matUG ( notfixeddofs, : ) = Ulib ;
fixednodes = nodalSprings(:,1) ;
gdlfixed = sort( nodes2dofs(fixednodes, ndofpnode) ) ;
Reactions = ((KG-KS)*matUG - matFs)(gdlfixed,:) ;

matFs ( fixeddofs, : ) = Fcon ; 
matFs ( notfixeddofs, : ) = Flib ;

tSystemResolution = toc ;
% ==============================================================================

% ==============================================================================
% -----------------------------    Post process   ------------------------------
tic ;
%~ nTimeSteps = size(matUG,2) ;

% Beams
ElemsSolic      = zeros(nelems,2*ndofpnode, nTimeSteps ) ;
dispsElemsMat   = zeros(nelems,2*ndofpnode, nTimeSteps ) ;

% Truss
trussDisps      = zeros(ntruss,2, nTimeSteps ) ;
trussStrain     = zeros(ntruss,1, nTimeSteps ) ;
normalForce     = zeros(ntruss+nbeam,1, nTimeSteps ) ;

% Plate
dispsPlateMat 	= zeros(nplate, 4*3, nTimeSteps) ;

platedofs = [ 2 4 5 ] ;

LocAxialdofs  = [ 1  1+ndofpnode                        ] ;
LocTorsndofs  = [ 2  2+ndofpnode                        ] ;
LocBendXYdofs = [ 3  6           3+ndofpnode 6+ndofpnode] ;
LocBendXZdofs = [ 5  4           5+ndofpnode 4+ndofpnode] ;

% Tetrahedron
strains  = [] ;
stresses = [] ;
cellDisp   = cell ;
% ----------------------------------------------------------------------




% Definition of matUts
if length(userLoadsFilename) > 0
  matUts = matUG ; 
else
  if norm(variableFext) > 0
    for currTime = 1:nTimeSteps
      currLoadFactor = loadFactorsFunc(currTime) ;
      if norm(constantFext) == 0
        matUts = [ matUts (matUG(:,1)*currLoadFactor) ] ;
      else
        matUts = [ matUts (matUG(:,1) + matUG(:,2)*currLoadFactor) ] ;
      end  
      %~ timeIndex = timeIndex +1 ;
    end
  else
    matUts = matUG ;
  end
end

% dispsElemsMat
for i = 1:nelems
  nodeselem = Conec(i,1:2)' ;
  dofselem  = nodes2dofs( nodeselem , ndofpnode ) ;
  for j = 1:nTimeSteps
    dispsElemsMat( i, :, j ) = matUts( dofselem, j )' ;
  end
end

% Postprocess for displacements and solicitations

for currTime = 1:nTimeSteps
	if currTime > 2
		if (size(matUts,2) == 1), error('numericalMethodParams vector must not be defined.'), end ;
  end
  for i = 1:nelems
    m = indexesElems(i) ;
    
    if Conec(i,7) == 1
    
      nodeselem  = Conec(i,1:2)' ;
      dofselem   = nodes2dofs( nodeselem , ndofpnode ) ;
      globalUelem = matUts(dofselem, currTime) ;
      l = ElemLengths(m) ;
      A = secGeomProps ( Conec(i,6), 1 ) ;
      E  = hyperElasParams{ Conec(i,5)}(1+1) ; 
      
      R = RotationMatrix ( ndofpnode, Local2GlobalMats{m} ) ;
      localUelem = R' * globalUelem ;
      
      trussDisps(m,:,currTime) = [localUelem(1) localUelem(1+6)] ; 
      trussStrain(m,currTime) = [ (trussDisps(m,2,currTime)-trussDisps(m,1,currTime))/l ] ;
      normalForce(m,currTime) = trussStrain(m,currTime) * E * A ;
      
    elseif Conec(i,7) == 2
    % obtains nodes and dofs of element
      nodeselem  = Conec(i,1:2)' ;
      dofselem   = nodes2dofs( nodeselem , ndofpnode ) ;
      globalUelem = matUts(dofselem, currTime) ; 
      l = ElemLengths(m) ;
      R = RotationMatrix ( ndofpnode, Local2GlobalMats{m} ) ;

      localUelem = R' * globalUelem ;


      RXYXZ = eye(4) ; RXYXZ(2,2) = -1; RXYXZ(4,4) = -1;   
      % Bending My - plane XZ
      if     elemReleases(m,1) == 1 && elemReleases(m,2) == 0
        if currTime == 1
          localUelem( 4) = -1* [  -3/(2*l)   0     3/(2*l) -0.5 ] * RXYXZ * localUelem(LocBendXZdofs) + elemSelfWAux(m,4)/((E*Iy/l^3)*4*l^2) ;
        else
          localUelem( 4) = -1* [  -3/(2*l)   0     3/(2*l) -0.5 ] * RXYXZ * localUelem(LocBendXZdofs) + elemUnifAux(m,4)/((E*Iy/l^3)*4*l^2) ;
        end   
      elseif elemReleases(m,1) == 0 && elemReleases(m,2) == 1
        if currTime == 1
          localUelem(10) = -1* [  -3/(2*l)  -0.5   3/(2*l)  0   ] * RXYXZ * localUelem(LocBendXZdofs) + elemSelfWAux(m,10)/((E*Iy/l^3)*4*l^2) ; 
        else
          localUelem(10) = -1* [  -3/(2*l)  -0.5   3/(2*l)  0   ] * RXYXZ * localUelem(LocBendXZdofs) + elemUnifAux(m,10)/((E*Iy/l^3)*4*l^2) ;
        end
      end
      % Bending Mz - plane XY
      if     elemReleases(m,3) == 1 && elemReleases(m,4) == 0
        if currTime == 1
          localUelem( 6) = +1 * [  -3/(2*l)   0     3/(2*l) -0.5 ] * localUelem(LocBendXYdofs) + elemSelfWAux(m,6)/((E*Iz/l^3)*4*l^2) ; 
        else
          localUelem( 6) = +1 * [  -3/(2*l)   0     3/(2*l) -0.5 ] * localUelem(LocBendXYdofs) + elemUnifAux(m,6)/((E*Iz/l^3)*4*l^2) ;
        end  
      elseif elemReleases(m,3) == 0 && elemReleases(m,4) == 1
        if currTime == 1
          localUelem(12) = +1 * [  -3/(2*l)  -0.5   3/(2*l)  0   ] * localUelem(LocBendXYdofs) + elemSelfWAux(m,12)/((E*Iz/l^3)*4*l^2) ; 
        else
          localUelem(12) = +1 * [  -3/(2*l)  -0.5   3/(2*l)  0   ] * localUelem(LocBendXYdofs) + elemUnifAux(m,12)/((E*Iz/l^3)*4*l^2) ; 
        end
      end
      dispsElemsMat( i, :, currTime ) = R * localUelem ;
      
      

    % -------------------------
    % Solicitations calculations
      KGelem = ElemKGs{m} ;
            
      if currTime == 1
				if selfWeightBoolean
					FGelem = KGelem * dispsElemsMat(i,:,currTime)' - elemSelfWLoc(i,:)' ;
				else 
					FGelem = KGelem * dispsElemsMat(i,:,currTime)' - elemUnifLoc(i,:)' ;
				end
			else
				FGelem = KGelem * dispsElemsMat(i,:,currTime)' - elemUnifLoc(i,:)' ;
      end
      
      %~ FGelem = KGelem * dispsElemsMat( i, :, currTime )'   ;
      FLelem = R' * FGelem ;
      
      ElemsSolic(i, :, currTime) = FLelem' ;
      
      epsxstrain = ( localUelem(1+6)-localUelem(1) ) / l ;
      normalForce(m,currTime) = epsxstrain * E * A ; 
    % -------------------------
    
    elseif Conec(i,7) == 3
    
      nodeselem =  Conec( i ,1:4) ;
      dofselem = nodes2dofs( nodeselem , ndofpnode ) ;
      dofstet = dofselem(1:2:end) ;
    
      strains ( m , : ) =            ( BMat{m} * matUts( dofstet, currTime ) )' ;
      stresses( m , : ) =  ( ConsMat{m} * strains(m,:)' )' ;

      cellDisp{1}   = 1 ;
      cellDisp{2}   = 'Displacements';
      cellDisp{3}   = matUts(1:2:end, currTime);

      cellStress{1} = 2 ;
      cellStress{2} = 'Stress';
      cellStress{3} = stresses ;
    elseif Conec(i,7) == 4 % plate element
      
      nodeselem = Conec ( i, 1:4 ) ;
      elemdofs = nodes2dofs( nodeselem, ndofpnode ) ;
			platedofs1 = [] ; % displacements dofs
      for j = 1:4
        if ismember(nodeselem(j), beamNodes)
          nodedofs = nodes2dofs(nodeselem(j), ndofpnode) ;
          platedofs1 = [ platedofs1 ; nodedofs([2 4 5]) ] ; % to fix
        else
          nodedofs = nodes2dofs(nodeselem(j), ndofpnode) ;
          platedofs1 = [ platedofs1 ; nodedofs([2 4 5]) ] ;
        end  
      end
       
      dispsPlateMat(i, :, currTime) = matUts(platedofs1, currTime) ;
      
    end
  end %endfor elements
  matNts = [ matNts normalForce(:,currTime) ] ;
end % endfor time

tSolicDisps = toc ;

% Analytic sol flag check

if analyticSolFlag ~= 0
  if analyticSolFlag == 2
    controlDisps = 0:deltaT:tf ;
    loadFactors = Reactions(analyticSolDofs,:);
    timesVec = 0:deltaT:tf ;
  elseif analyticSolFlag == 3
    controlDisps = matUts( analyticSolDofs, indexTimeSol ) ;
  elseif analyticSolFlag == 4   
    controlDisps = dispsElemsMat ;
  elseif analyticSolFlag == 5 
    controlDisps = Reactions ( analyticSolDofs, indexTimeSol ) ;
  end
else
	controlDisps = 0 ;
	loadFactors = 0 ;
end

%~ lw  = 2   ; ms  = 5.5 ;
%~ lw2 = 3.2 ; ms2 = 23 ;
%~ plotfontsize = 22 ;
%~ % Plot moment 
%~ hold on
%~ for i=1:nelems
    %~ plot3( coordsElemsMat(i,[1 7]), coordsElemsMat(i,[3 9]), coordsElemsMat(i,[5 11]),'b--o','linewidth',lw*0.8,'markersize',ms*0.8);
%~ end
%~ for i =1:nelems
%~ i=2


	%~ xsref = coordsElemsMat(i, [ 1   7  ] ) ;
	%~ ysref = coordsElemsMat(i, [ 1+2 7+2] ) ;
	%~ zsref = coordsElemsMat(i, [ 1+4 7+4] ) ;

	%~ nPlotPoints = 30 ;

	%~ l      = sqrt( sum( (xsref(2)-xsref(1))^2 + (ysref(2)-ysref(1))^2  + (zsref(2)-zsref(1))^2 ) ) ;

	%~ R = RotationMatrix ( ndofpnode, Local2GlobalMats{i} ) ;
	%~ localUelem = R' * dispsElemsMat(i,:,:)' ;

	%~ E  = hyperElasParams{ Conec(i,5)}(1+1) ; 
	%~ Iy = secGeomProps ( Conec(i,6), 2 ) ;  

	%~ xsloc  = linspace( 0 , l, nPlotPoints )' ;
	
	%~ M     = zeros( size(xsloc) ) ;
	%~ for j = 1:nPlotPoints
		%~ Rb = eye(4);
    %~ Rb(2,2) = -1;
    %~ Rb(4,4) = -1;
		%~ N = bendingInterFuns( xsloc(j), l, 2 );
		%~ M(j) = -N  * Rb' * (localUelem([ 5  4           5+ndofpnode 4+ndofpnode]) * 1* 1 ) + elemUnifLoc(i,[ 4])' + elemSelfWLoc(i, [4])';
	%~ end
	
		%~ exL = Local2GlobalMats{i}(:,1) ;
	%~ xsRef = ( xsref(1) + xsloc )' ;
	%~ if zsref(2) < zsref(1)
		%~ zsRef = ( zsref(1) - xsloc )' ;
	%~ else 
		%~ zsRef = ( zsref(1) + xsloc )' ;
	%~ end
	%~ MRef 	= (	zsref(1) + M ) ;
	%~ if isequal(abs(exL), [0 0 1]') 
		%~ i
		%~ plot3( xsref(1)+M, zeros(size(xsloc)), zsRef, 'b-', 'linewidth', lw, 'markersize', ms );
	%~ else
		%~ plot3( xsRef, zeros(size(xsloc)), zsref(1)+M, 'b-', 'linewidth', lw, 'markersize', ms );
	%~ end
	%~ grid on
	%~ view( [1 -1 1])
%~ end

% ==============================================================================
if plotParamsVector(1)>0
  fprintf(' done.\n');
end
