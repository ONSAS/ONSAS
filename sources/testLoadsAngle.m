clear all, close all

%~ indselim = unique( round( rand(round(0.5*n),1)*(n*.95) + 1 ) )
%~ ts (indselim) = [];

%~ Rx = 4;
%~ Ry = 2;
%~ xs = linspace(0,1,n)';
%~ ys = xs.^2 ;
%~ xs = Rx*cos(ts);
%~ ys = Ry*sin(ts);
%~ zs = ts*0;


%~ EI  =1 ;

%~ Nodes = [ xs ys zeros(size(ys)) ];
Nodes = [ 0 0 0 ;
          1 1 0 ;
          2 0 0 ];

Conec = [ 1 2; 2 3 ];

bendStiff = [0 1 0 ]

fextAngSpr = loadsAngleSpring( Nodes, Conec, bendStiff )



%~ 1/ ( Rx^2 / Ry )
%~ 1/ ( Ry^2 / Rx )


%~ figure
%~ plot(xs,ys,'o-g')
%~ plot3(xs,ys,zs,'o-g')
%~ hold on
%~ quiver(xs,ys,zs,fextAngSpr(1:6:end),fextAngSpr(3:6:end),fextAngSpr(3:6:end),'o-r')
%~ quiver3(xs,ys,zs,fextAngSpr(1:6:end),fextAngSpr(3:6:end),fextAngSpr(3:6:end),'o-r')
%~ axis equal

%~ figure
%~ plot(ts,chis,'x-r')
%~ Nodes



