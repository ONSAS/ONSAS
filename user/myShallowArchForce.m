function f = myShallowArchForce(t)

%Force and direction
nnodes  = 25 ;
nnodes2 =  22; %nodo a -L/3

Fzmax = 6 * 1e6    ;
OmegaFz = 10       ;

f = zeros( nnodes*6,1);

f( (nnodes2-1)*6 + 5 )  = Fzmax*sin(OmegaFz*t) ;
