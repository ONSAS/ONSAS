% Copyright (C) 2019, Jorge M. Perez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquin Viera, Mauricio Vanzulli  
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


function [nodalConstUnif, nodalSelfW, nodalPlate, elemUnifLoc, elemSelfWLoc, elemUnifAux, elemSelfWAux] = assemblyUniform(unifLoad, selfWeightBoolean, hyperElasParams, Conec, Nodes, indexesElems, ElemLengths, Local2GlobalMats, elemReleases, secGeomProps, Lx, Ly)  

	ndofpnode = 6 ;
	nelems = size(Conec,1) ;
	nnodes = size(Nodes,1) ;
	nodalConstUnif = zeros(ndofpnode*nnodes,1) ;
	nodalSelfW = zeros(ndofpnode*nnodes,1) ;
	nodalPlate = zeros(ndofpnode*nnodes,1) ;
	% Para sumar a los releases, carga distribuida comun
	elemUnifAux = zeros(nelems,ndofpnode*2) ;
	elemSelfWAux = zeros(nelems,ndofpnode*2) ;
	% Para restar a FGelem, carga distribuida modificado segun release
	elemUnifLoc = zeros(nelems, ndofpnode*2) ;
	elemSelfWLoc = zeros(nelems, ndofpnode*2) ;

	for i = 1:nelems
		m = indexesElems(i) ;
		
		if Conec(i,7) == 1
		
			if selfWeightBoolean == 1
				A = secGeomProps ( Conec(i,6), 1 ) ;
				l = ElemLengths(m) ;
        rho = hyperElasParams{Conec(i,5)}(end) ;
				R = RotationMatrix ( ndofpnode, Local2GlobalMats{m} ) ;
				exL = Local2GlobalMats{m}(:,1); eyL = Local2GlobalMats{m}(:,2); ezL = Local2GlobalMats{m}(:,3) ;
				nodeselem = Conec(i,1:2)' ;
				aux = nodes2dofs ( nodeselem, ndofpnode ) ;
				
				cosZZ = ezL' * [0 0 1]'; senXZ = exL'*[0 0 1]' ;
				
				qW = rho * A ;
				qlloc = qW * l * -[ 1/2  0 0   0    0  0 1/2 0 0   0   0  0 ]' ;
				qploc = qW * l * -[ 	0	 0 0	 0	1/2	 0	 0 0 0 	 0 1/2	0 ]' ;
				ql = qlloc * senXZ ;
				qp = qploc * cosZZ ;
				qpAux = qp ;
				
				elemSelfWAux(i,:) = (R * ( sum(qpAux,2) + sum(ql,2) ))' + elemSelfWAux(i,:) ;
				elemSelfWLoc(i,:) = (R * ( sum(ql,2) 		+ sum(qp,2) ))' + elemSelfWLoc(i,:) ; 
				nodalSelfW(aux) = elemSelfWLoc(i,:)' + nodalSelfW(aux) ;
				
				
			end
		
		elseif Conec(i,7) == 2
			A = secGeomProps ( Conec(i,6), 1 ) ;
			l = ElemLengths(m) ;
			R = RotationMatrix ( ndofpnode, Local2GlobalMats{m} ) ;
			exL = Local2GlobalMats{m}(:,1); eyL = Local2GlobalMats{m}(:,2); ezL = Local2GlobalMats{m}(:,3) ;
			nodeselem = Conec(i,1:2)' ;
			aux = nodes2dofs ( nodeselem, ndofpnode ) ;
			
			if ismember(i,unifLoad)
					
				nPosElem = find(unifLoad(:,1) == i) ;
				vecqx = unifLoad(nPosElem, 3); vecqy = unifLoad(nPosElem, 4); vecqz = unifLoad(nPosElem, 5) ;
				
				for j = 1:length(nPosElem)
					qx = vecqx(j); qy = vecqy(j); qz = vecqz(j) ;
					
					% local uniform forces
					qlloc =  l * [ 1/2  0 0   0    0  0 1/2 0 0   0   0  0 ]' ;
					qpAuxloc = l * [ 0 0 		0 	-l/12 	1/2 		0 	0 	0 		0 	l/12 	1/2 		0 ; ...
													 0 0 	1/2 			0 		0	 l/12 	0 	0 	1/2 	 	 0 		0 -l/12 ; ...
													 0 0 		0 	-l/12 	1/2 		0 	0 	0 		0 	l/12 	1/2 		0 ]' ;
					if elemReleases(m,1) == 1 && elemReleases(m,2) == 0
						qploc =  l * [ 0 0 		0 			0 	3/8 		0  	0 	0 		0 	 l/8 	5/8 		0 ; ...
													 0 0 	3/8 			0 		0 		0  	0 	0 	5/8 		 0 		0  -l/8 ; ...
													 0 0 		0 			0 	3/8 		0  	0 	0 		0 	 l/8 	5/8 		0 ]' ;					   
					elseif elemReleases(m,1) == 0 && elemReleases(m,2) == 1
						qploc =  l * [ 0 0 		0  	 -l/8 	5/8 		0  	0 	0 		0 		 0 	3/8 		0 ; ...
													 0 0 	5/8 			0 		0 	l/8  	0 	0 	3/8 		 0 		0 		0 ; ...
													 0 0 		0  	 -l/8 	5/8 		0  	0 	0 		0 		 0 	3/8 		0 ]' ;
					elseif elemReleases(m,1) == 0 && elemReleases(m,2) == 0
						qploc =  l * [ 0 0 		0 	-l/12 	1/2 		0 	0 	0 		0 	l/12 	1/2 		0 ; ...
													 0 0 	1/2 			0 		0	 l/12 	0 	0 	1/2 	 	 0 		0 -l/12 ; ...
													 0 0 		0 	-l/12 	1/2 		0 	0 	0 		0 	l/12 	1/2 		0 ]' ; 
					end
					
					if unifLoad(nPosElem(j),2) == 1 % GLOBAL AXES
						% X global
						cosZX = ezL' * [1 0 0]'; senXZ = exL'*[0 0 1]' ;
						% Y global
						cosYY = eyL' * [0 1 0]'; senXY = exL'*[0 1 0]' ;
						% Z global
						cosZZ = ezL' * [0 0 1]'; senXZ = exL'*[0 0 1]' ;
						
						ql = qlloc * diag([qx qy qz]' * [senXZ senXY senXZ])' ;
						qp = qploc * diag([qx qy qz]' * [cosZX cosYY cosZZ]) ;
						qpAux = qpAuxloc * diag([qx qy qz]' * [cosZX cosYY cosZZ]) ;
						
					elseif unifLoad(nPosElem(j),2) == 0 % LOCAL AXES
						
						ql = qx * qlloc ;
						qp = qploc * [0 qy qz]' ;
						qpAux = qpAuxloc * [0 qy qz]' ;
						
					end % endif GLOBAL/LOCAL
					elemUnifAux(i,:) = (R * ( sum(qpAux,2) 	+ sum(ql,2) ))' + elemUnifAux(i,:) ;
					elemUnifLoc(i,:) = (R * ( sum(ql,2) 		+ sum(qp,2) ))' + elemUnifLoc(i,:) ;
					
					nodalConstUnif(aux) = elemUnifLoc(i,:)' + nodalConstUnif(aux) ;
				end % endfor nposElem
			
      end % endif ismember unifLoad
      
      if selfWeightBoolean
				rho = hyperElasParams{Conec(i,5)}(end) ;
				qW = rho * A ;
				cosZZ = ezL' * [0 0 1]'; senXZ = exL'*[0 0 1]' ;
				qlloc = qW * l * -[ 1/2  0 0   0    0  0 1/2 0 0   0   0  0 ]' ;
				qpAuxloc = qW * l * -[ 0 0 0 -l/12 	1/2 	0 0 0 0  l/12 1/2 0 ]' ;
				if elemReleases(m,1) == 1 && elemReleases(m,2) == 0
          qploc = qW * l * -[ 0 0 0 	  0 	3/8 	0 0 0 0  	l/8 5/8 0 ]' ;
        elseif elemReleases(m,1) == 0 && elemReleases(m,2) == 1
          qploc = qW * l * -[ 0 0 0  -l/8 	5/8 	0 0 0 0 	 	0 3/8 0 ]' ;
        elseif elemReleases(m,1) == 0 && elemReleases(m,2) == 0
					qploc = qW * l * -[ 0 0 0 -l/12 	1/2 	0 0 0 0  l/12 1/2 0 ]' ;
        end
				ql = qlloc * senXZ ;
				qp = qploc * cosZZ ;
				qpAux = qpAuxloc * cosZZ ;
				
				elemSelfWAux(i,:) = (R * ( sum(qpAux,2) + sum(ql,2) ))' + elemSelfWAux(i,:) ;
				elemSelfWLoc(i,:) = (R * ( sum(ql,2) 		+ sum(qp,2) ))' + elemSelfWLoc(i,:) ; 
				nodalSelfW(aux) = elemSelfWLoc(i,:)' + nodalSelfW(aux) ;
				
      end 
      
		elseif Conec(i,7) == 4
			nodeselem = Conec(i,1:4) ;
			elemdofs = nodes2dofs(nodeselem,ndofpnode) ;
			
			a = Lx(m) / 2 ;
			b = Ly(m) / 2 ;
			q = unifLoad(i,5) ;

			qelem = 4*q*a*b* [ 0 a/12 0 b/12 1/4 0 0 -a/12 0 b/12 1/4 0 0 -a/12 0 -b/12 1/4 0 0 a/12 0 -b/12 1/4 0 ]' ;
      
			nodalPlate(elemdofs) = nodalPlate(elemdofs) + qelem ;
      
      if selfWeightBoolean
        rho = hyperElasParams{Conec(i,5)}(end) ;
        % todo
      end  
        
		end
	end
end
