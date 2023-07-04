
% implementation of the DKT plate triangle element based on https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620210709

function [ fs, ks, stress, strain ] = internal_forces_plate_triangle( ...
  elemCoords, elemDisps, hyperElasModel, elemConstitutiveParams, paramOut, t, planeStateFlag, dotdotdispsElem, density, previous_state )

E = 200e9;
nu = .3;
t = .01 ;

elemCoords = [0 0 0 1 .5 0 0 2 0];

  assert( norm( elemCoords(3:3:end))==0,'only xy plates are considered' )

  x1 = elemCoords(1);
  y1 = elemCoords(2);

  x2 = elemCoords(4);
  y2 = elemCoords(5);

  x3 = elemCoords(7);
  y3 = elemCoords(8);

mat_cross = [ x2-x3 x3-x1; y2-y3 y3-y1; 0 0 ]

norm( cross( mat_cross(:,1), mat_cross(:,2)) )

det( mat_cross(1:2,:) )

B=[y2-y3; y3-y1; y1-y2];
C=[x3-x2; x1-x3; x2-x1];

DET = ( B(1) * C(2) - B(2) * C(1) ) * 24

% isotropic
D = t^3/12 * E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)*.5 ] ;

PP = [ 12 4 4 ; 4 2 1; 4 1 2 ];

DD =  zeros(9,9);

i = 1; j = 1;  DD(nodes2dofs(i,3),nodes2dofs(j,3)) = D(i,j)*PP/DET ;
i = 1; j = 2;  DD(nodes2dofs(i,3),nodes2dofs(j,3)) = D(i,j)*PP/DET ;
i = 2; j = 2;  DD(nodes2dofs(i,3),nodes2dofs(j,3)) = D(i,j)*PP/DET ;
i = 3; j = 3;  DD(nodes2dofs(i,3),nodes2dofs(j,3)) = D(i,j)*PP/DET ;

i = 2; j = 1;  DD(nodes2dofs(i,3),nodes2dofs(j,3)) = D(i,j)*PP'/DET ;


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

QQ = zeros(9,9);

QQ(1,:) = B(2)*GG(1,:) + B(3)*GG(2,:);
QQ(2,:) = 2*B(2)*GG(3,:) + B(3)*GG(4,:);
QQ(3,:) = B(2)*GG(4,:) + 2*B(3)*GG(5,:);

QQ(4,:) = -C(2)*GG(6,:) - C(3)*GG(7,:);
QQ(5,:) = -2*C(2)*GG(8,:) - C(3)*GG(9,:) ;
QQ(6,:) = -C(2)*GG(9,:) -2*C(3)*GG(10,:) ;

QQ(7,:) =  C(2)*GG(1,:) + C(3)*GG(2,:) ...
          -B(2)*GG(6,:) - B(3)*GG(7,:) ;

QQ(8,:) = 2*C(2)*GG(3,:) + C(3)*GG(4,:) ...
         -2*B(2)*GG(8,:) - B(3)*GG(9,:) ;

QQ(9,:) = C(2)*GG(4,:) + 2*C(3)*GG(5,:) ...
         -B(2)*GG(9,:) - 2*B(3)*GG(10,:) ;

K = QQ' * DD * QQ ;


Kred = K(7:9,7:9) ;

u = Kred \ [ -1e3 ; 0 ; 0]