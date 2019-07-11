clear all, close all

%~ function fintAngSpr
n=40;

ts = linspace(0,pi,n)' ;

indselim = unique( round( rand(round(0.5*n),1)*(n*.95) + 1 ) )
%~ ts (indselim) = [];


Rx = 4;
Ry = 2;
%~ xs = linspace(0,1,n)';
%~ ys = xs.^2 ;
xs = Rx*cos(ts);
ys = Ry*sin(ts);
zs = ts*0;


EI  =1 ;

Nodes = [ xs ys zeros(size(ys)) ];

n=size(Nodes,1);

fextAngSpr = zeros( 6*n,1);

          
Conec = [ (1:(n-1))' (2:(n))' ];

chis = zeros( size(Nodes,1) ,1 );          

for i=1 : (size(Conec,1) -1 )
  indxim1 = Conec(i  ,1) ;
  indxi   = Conec(i  ,2) ;
  indxip1 = Conec(i+1,2) ;
  
  ti   = Nodes(indxi  ,:) - Nodes(indxim1,:) ;
  tip1 = Nodes(indxip1,:) - Nodes(indxi  ,:) ;

  li   = norm( ti);
  lip1 = norm( tip1);
  
  dtds = ( tip1/lip1 - ti/li ) / ( lip1/2 + li/2 ) ;
  chis(i+1) = norm( dtds );
  %~ normal = dtds / chis(i+1) ;

  Fiim1 = -dtds *  EI * chis(i+1) / norm( cross( ti  , dtds) )  ;
  Fim1i =  dtds *  EI * chis(i+1) / norm( cross( ti  , dtds) )  ;
  
  Fiip1 = -dtds *  EI * chis(i+1) / norm( cross( tip1, dtds) )  ;
  Fip1i =  dtds *  EI * chis(i+1) / norm( cross( tip1, dtds) )  ;

  aux = nodes2dofs( indxim1,6) ;
  fextAngSpr( aux(1:2:6) ) = Fim1i ; 

  aux = nodes2dofs( indxi  ,6) ;
  fextAngSpr( aux(1:2:6) ) = Fiim1 + Fiip1 ; 

  aux = nodes2dofs( indxip1,6) ;
  fextAngSpr( aux(1:2:6) ) = Fip1i ; 
  
  %~ pause
end



1/ ( Rx^2 / Ry )
1/ ( Ry^2 / Rx )


figure
%~ plot(xs,ys,'o-g')
plot3(xs,ys,zs,'o-g')
hold on
%~ quiver(xs,ys,zs,fextAngSpr(1:6:end),fextAngSpr(3:6:end),fextAngSpr(3:6:end),'o-r')
quiver3(xs,ys,zs,fextAngSpr(1:6:end),fextAngSpr(3:6:end),fextAngSpr(3:6:end),'o-r')
axis equal

figure
plot(ts,chis,'x-r')
%~ Nodes



