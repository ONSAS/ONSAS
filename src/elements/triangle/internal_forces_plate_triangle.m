
% implementation of the DKT plate triangle element based on https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620210709

function [ fs, ks ] = internal_forces_plate_triangle( elemCoords, elemDisps, modelName, modelParams, thickness )

% assertions
assert( norm( elemCoords(3:3:end))==0, 'only plates in x-y plane are considered' )
assert( strcmp( modelName, 'elastic-linear'), ' linear elastic model is implemented' )

E = modelParams(1)  ;  nu = modelParams(2) ;

x1G = elemCoords(0*3+1);  y1G = elemCoords(0*3+2);
x2G = elemCoords(1*3+1);  y2G = elemCoords(1*3+2);
x3G = elemCoords(2*3+1);  y3G = elemCoords(2*3+2);

x1 = x1G ; y1 = y1G ; 
x2 = x2G ; y2 = y2G ; 
x3 = x3G ; y3 = y3G ; 

mat_cross = [ x2-x3 x3-x1; y2-y3 y3-y1; 0 0 ] ;

det( mat_cross(1:2,:) );

B=[y2-y3; y3-y1; y1-y2];
C=[x3-x2; x1-x3; x2-x1];

DET = ( B(1) * C(2) - B(2) * C(1) ) * 24 ;

% isotropic
D = thickness^3/12 * E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)*.5 ] ;

PP = 1/24 * [ 12 4 4 ; 4 2 1; 4 1 2 ];
  
DD =  zeros(9,9);
i = 1; j = 1;  DD( nodes2dofs(i,3), nodes2dofs(j,3)) = D(i,j)*PP/DET ;
i = 1; j = 2;  DD( nodes2dofs(i,3), nodes2dofs(j,3)) = D(i,j)*PP/DET ;
i = 2; j = 2;  DD( nodes2dofs(i,3), nodes2dofs(j,3)) = D(i,j)*PP/DET ;
i = 3; j = 3;  DD( nodes2dofs(i,3), nodes2dofs(j,3)) = D(i,j)*PP/DET ;
i = 2; j = 1;  DD( nodes2dofs(i,3), nodes2dofs(j,3)) = D(i,j)*PP/DET ;
  
ALS = B.^2 + C.^2;
PT = [ (6*C./ALS)' ; (6*B./ALS)'] ;
RS = [ (3*(C.^2)./ALS)' ; (3*(B.^2)./ALS)' ];
Q = 3 * B .* C ./ ALS ;

GG = zeros(10,9) ;
KOD = [ 1  1 ; 2 3; 3 2; 4 4; 5 6; 6 5; 7 7; 8 9 ; 9 8]' ;

for I=1:2
  II=(I-1)*5 ;
    
  GG(II+1, KOD(I,1)) =  PT(I, 3);
  GG(II+2, KOD(I,1)) = -PT(I, 2);
  GG(II+3, KOD(I,1)) = -PT(I, 3);
  GG(II+4, KOD(I,1)) =  PT(I, 2)-PT(I,3);
  GG(II+5, KOD(I,1)) =  PT(I, 2);
  
  GG(II+1, KOD(I,2)) = -Q(3);
  GG(II+2, KOD(I,2)) = -Q(2);
  GG(II+3, KOD(I,2)) =  Q(3);
  GG(II+4, KOD(I,2)) =  Q(2)+Q(3);
  GG(II+5, KOD(I,2)) =  Q(2);

  GG(II+1, KOD(I,3)) = -1 -RS(I,3);
  GG(II+2, KOD(I,3)) = -1 -RS(I,2);
  GG(II+3, KOD(I,3)) = RS(I,3);
  GG(II+4, KOD(I,3)) = RS(I,2)+RS(I,3);
  GG(II+5, KOD(I,3)) = RS(I,2);

  GG(II+1, KOD(I,4)) = -PT(I,3);
  GG(II+3, KOD(I,4)) =  PT(I,3);
  GG(II+4, KOD(I,4)) =  PT(I,1)+PT( I ,3);
  
  GG(II+1, KOD(I,5)) = -Q(3);
  GG(II+3, KOD(I,5)) =  Q(3);
  GG(II+4, KOD(I,5)) =  Q(3)-Q(1);
  
  GG(II+1, KOD(I,6)) = 1-RS(I,3)       ;
  GG(II+3, KOD(I,6)) = RS(I,3)         ;
  GG(II+4, KOD(I,6)) = RS(I,3) -RS(I,1);
  
  GG(II+2, KOD(I,7)) =  PT(I,2)        ;
  GG(II+4, KOD(I,7)) = -PT(I,1)-PT(I,2);
  GG(II+5, KOD(I,7)) = -PT(I,2)        ;
  
  GG(II+2, KOD(I,8)) = -Q(2)           ;
  GG(II+4, KOD(I,8)) =  Q(2)-Q(1)      ;
  GG(II+5, KOD(I,8)) =  Q(2)           ;
  
  GG(II+2, KOD(I,9)) = 1-RS(I,2)       ;
  GG(II+4, KOD(I,9)) = RS(I,2)-RS(I,1) ;
  GG(II+5, KOD(I,9)) = RS(I,2)         ;
end

matX = [ B(2)*GG(1,:)+B(3)*GG(2,:)    ;
       2*B(2)*GG(3,:)+B(3)*GG(4,:)    ;
         B(2)*GG(4,:)+2*B(3)*GG(5,:)  ] ;
matY = [-C(2)*GG(6,:)-C(3)*GG(7,:)    ;
        -2*C(2)*GG(8,:)-C(3)*GG(9,:)  ;
        -C(2)*GG(9,:)-2*C(3)*GG(10,:) ] ;
matZ = [ C(2)*GG(1,:)+C(3)*GG(2,:)-B(2)*GG(6,:)-B(3)*GG(7,:)        ;
         2*C(2)*GG(3,:)+C(3)*GG(4,:)-2*B(2)*GG(8,:)-B(3)*GG(9,:)    ;
         C(2)*GG(4,:)+2*C(3)*GG(5,:)-B(2)*GG(9,:)-2*B(3)*GG(10,:) ] ;

XYZ = [ matX; matY; matZ ];

K = XYZ' * DD * XYZ ;

gamma_vec = [ 1 1 0 ] ;
aux = zeros(3, 9) ;
aux(1,0+(1:3)) = gamma_vec ;
aux(2,3+(1:3)) = gamma_vec ;
aux(3,6+(1:3)) = gamma_vec ;
vec_xis = 1/DET* aux * ( XYZ * elemDisps ) ;
Moments = D * vec_xis ;
Mx = Moments(1) 

fint = K * elemDisps ;

fs = {fint};  ks = {K};

% fout = zeros(18,1);
% fout(dofs_sort) = fint ;

% Kout = zeros(18,18) ;
% Kout(dofs_sort,dofs_sort) = K ;

% fs = {fout};  ks = {Kout};