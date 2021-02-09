function [NodesCoord, Conec] =  createSolidRegularMesh( ndivs, Lx, Ly, Lz )
nx = ndivs(1);
ny = ndivs(2);
nz = ndivs(3);

nnodes = (nx+1)*(ny+1)*(nz+1);
nelems = prod( ndivs)*6;

NodesCoord = zeros( nnodes, 3 ) ;

coodrsx = linspace(0,Lx, nx+1)';
coodrsy = linspace(0,Ly, ny+1)';
coodrsz = linspace(0,Lz, nz+1)';

for iz=1:(nz+1)
  for iy=1:(ny+1)
    prevNodes = (iy-1)*(nx+1) + (iz-1)*(nx+1)*(ny+1)
    NodesCoord( prevNodes + (1:(nx+1)),:) = [ coodrsx coodrsy(iy)*ones(nx+1,1) coodrsz(iz)*ones(nx+1,1) ] ;
  end
end

Conec = zeros( nelems, 4) ;
cont = 0;

for iz=1:nz
  for iy=1:ny
    for ix=1:nx
      n1 = ix + (iy-1)*(nx+1) + (iz-1)*(nx+1)*(ny+1) ;
      n2 = n1 + (nx+1)*(ny+1) ;
      n3 = n1 + (nx+1)*(ny+1) + (nx+1) ;
      n4 = n1 + (nx+1) ;
      
      n5 = n1 + 1 ;
      n6 = n2 + 1 ;
      n7 = n3 + 1 ;
      n8 = n4 + 1 ;
      
      Conec( cont + (1:6),:) = [ n1 n4 n2 n6 ; ...
                                 n6 n2 n3 n4 ; ... 
                                 n4 n3 n6 n7 ; ... 
                                 n4 n1 n5 n6 ; ... 
                                 n4 n6 n5 n8 ; ...
                                 n4 n7 n6 n8 ] ;
      cont = cont + 6 ;
    end
  end
end


