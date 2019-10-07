%function for assembly of tangent stiffness matrix and/or internal forces vector. Input parameter used to set output: only internal forces vector (1) or only tangent matrices (2)  

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


function [FintGt, KT, StrainVec, StressVec ] = assemblyFintVecTangMat ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Ut, bendStiff, paramOut )

ndofpnode = 6;
nelems = size(Conec,1);

KT  = sparse( length(Ut) , length(Ut)  ) ;
FintGt = zeros( length(Ut) ,1 ) ;

dispsElemsMat = zeros(nelems,2*ndofpnode) ;
for i=1:nelems
  % obtains nodes and dofs of element
  nodeselem = Conec(i,1:2)' ;
  dofselem  = nodes2dofs( nodeselem , ndofpnode ) ;
  dispsElemsMat( i, : ) = Ut(dofselem)' ;
end


StrainVec   = zeros( nelems, 6 ) ;
StressVec   = zeros( nelems, 6 ) ;

for elem = 1:nelems
  % obtains nodes and dofs of element
  nodeselem = Conec(elem,1:2)' ;
  dofselem  = nodes2dofs( nodeselem , ndofpnode ) ;

  switch Conec(elem,7)

  % -------------------------------------------
  case 1 % Co-rotational Truss

    sizeTensor = 1 ;

    A  = secGeomProps(Conec(elem,6),1) ;

    if paramOut == 2
      [ KTe, KL0e ] = elementTruss3DTangentMats( ...
        coordsElemsMat(elem,:), dispsElemsMat(elem,:), hyperElasParamsMat( Conec(elem,5),:) , A );
    else
      [ Finte, stress, dstressdeps, strain ] = elementTruss3DInternLoads( ...
        coordsElemsMat(elem,:), dispsElemsMat(elem,:), hyperElasParamsMat( Conec(elem,5),:) , A );
    end

  % -------------------------------------------
  case 2 % Co-rotational Frame element (bernoulli beam)
    sizeTensor = 1 ;

    A   = secGeomProps(Conec(elem,6),1) ;
    Iyy = secGeomProps(Conec(elem,6),2) ;
    Izz = secGeomProps(Conec(elem,6),3) ;
    J   = secGeomProps(Conec(elem,6),4) ;

    xs = coordsElemsMat(elem,1:2:end)'        ;
    E  = hyperElasParamsMat( Conec(elem,5),2) ;
    nu = hyperElasParamsMat( Conec(elem,5),3) ;
    G  = E/(2*(1+nu)) ;

    [ Finte, KTe, strain, stress ]= elementBeam3DInternLoads( xs, dispsElemsMat(elem,:)' , [E G A Iyy Izz J] ) ;

    KL0e = KTe;

  % -------------------------------------------
  case 3 % linear solid element
    
    sizeTensor = 6 ;

    [ Finte, KTe, strain, stress ]= elementBeam3DInternLoads( xs, dispsElemsMat(elem,:)' , [E G A Iyy Izz J] ) ;

  end


  % -------------------------------------------
  if paramOut == 1
    % internal loads vector assembly
    FintGt ( dofselem ) = FintGt( dofselem ) + Finte ;
  
  	StrainVec(elem,(1:sizeTensor) ) = strain ;
		StressVec(elem,(1:sizeTensor) ) = stress ;
  
  elseif paramOut == 2
    % matrices assembly
    KT  (dofselem,dofselem) = KT  (dofselem,dofselem) + KTe     ;
  end

end

KT     = KT  + KS ;
FintGt = FintGt + KS*Ut ;

if length(bendStiff) >0

  Nodes = conv ( Conec, coordsElemsMat+dispsElemsMat ) ;

  [ ~, KTAngSpr ] = loadsAngleSpring( Nodes, Conec, bendStiff ) ;

  fextAngSpr = KTAngSpr*Ut ;

  KT     += KTAngSpr   ;
  FintGt += fextAngSpr ;

end

% ------------------------------------



function nodesmat = conv ( conec, coordsElemsMat ) 

nodesmat  = [] ;
nodesread = [] ;

for i=1:size(conec,1)
  for j=1:2
    if length( find( nodesread == conec(i,j) ) ) == 0
      nodesmat( conec(i,j),:) = coordsElemsMat( i, (j-1)*6+(1:2:5) ) ;
    end
  end
end
