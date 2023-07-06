

function [ fs, ks ] = internal_forces_plate_triangle( elemCoords, elemDisps, hyperElasModel, hyperElasParams, thickness )

% assertions
assert( norm( elemCoords(3:3:end))==0, 'only xy plates are considered' )
assert( strcmp( hyperElasModel, 'linearElastic'), ' linear elastic model is implemented' )

E = hyperElasParams(1)  ;  nu = hyperElasParams(2) ;

x1 = elemCoords(1);  y1 = elemCoords(2);
x2 = elemCoords(4);  y2 = elemCoords(5);
x3 = elemCoords(7);  y3 = elemCoords(8);

mat_cross = [ x2-x3 x3-x1; y2-y3 y3-y1; 0 0 ] ;

% implementation of the DKT plate triangle element based on https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620210709

% norm( cross( mat_cross(:,1), mat_cross(:,2)) );

% det( mat_cross(1:2,:) );

% B=[y2-y3; y3-y1; y1-y2];
% C=[x3-x2; x1-x3; x2-x1];

% DET = ( B(1) * C(2) - B(2) * C(1) ) * 24 ;

% % isotropic
% D = thickness^3/12 * E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)*.5 ] ;

% PP = [ 12 4 4 ; 4 2 1; 4 1 2 ];

% DD =  zeros(9,9);
% i = 1; j = 1;  DD(nodes2dofs(i,3),nodes2dofs(j,3)) = D(i,j)*PP/DET ;
% i = 1; j = 2;  DD(nodes2dofs(i,3),nodes2dofs(j,3)) = D(i,j)*PP/DET ;
% i = 2; j = 2;  DD(nodes2dofs(i,3),nodes2dofs(j,3)) = D(i,j)*PP/DET ;
% i = 3; j = 3;  DD(nodes2dofs(i,3),nodes2dofs(j,3)) = D(i,j)*PP/DET ;
% i = 2; j = 1;  DD(nodes2dofs(i,3),nodes2dofs(j,3)) = D(i,j)*PP'/DET ;

% ALS = B.^2 + C.^2;
% PT = [ (6*C./ALS)' ; (6*B./ALS)'] ;
% RS = [ (3*(C.^2)./ALS)' ; (3*(B.^2)./ALS)' ];
% Q = 3 * B .* C ./ ALS ;

% GG = zeros(10,9) ;
% KOD = [ 1  1 ; 2 3; 3 2; 4 4; 5 6; 6 5; 7 7; 8 9 ; 9 8]' ;

% for I=1:2
%   II=(I-1)*5 ;
  
%   GG(II+1, KOD(I,1)) =  PT(I, 3);
%   GG(II+2, KOD(I,1)) = -PT(I, 2);
%   GG(II+3, KOD(I,1)) = -PT(I, 3);
%   GG(II+4, KOD(I,1)) =  PT(I, 2)-PT(I,3);
%   GG(II+5, KOD(I,1)) =  PT(I, 2);
  
%   GG(II+1, KOD(I,2)) = -Q(3);
%   GG(II+2, KOD(I,2)) = -Q(2);
%   GG(II+3, KOD(I,2)) =  Q(3);
%   GG(II+4, KOD(I,2)) =  Q(2)+Q(3);
%   GG(II+5, KOD(I,2)) =  Q(2);

%   GG(II+1, KOD(I,3)) = -1 -RS(I,3);
%   GG(II+2, KOD(I,3)) = -1 -RS(I,2);
%   GG(II+3, KOD(I,3)) = RS(I,3);
%   GG(II+4, KOD(I,3)) = RS(I,2)+RS(I,3);
%   GG(II+5, KOD(I,3)) = RS(I,2);

%   GG(II+1, KOD(I,4)) = -PT(I,3);
%   GG(II+3, KOD(I,4)) =  PT(I,3);
%   GG(II+4, KOD(I,4)) =  PT(I,1)+PT( I ,3);

%   GG(II+1, KOD(I,5)) = -Q(3);
%   GG(II+3, KOD(I,5)) =  Q(3);
%   GG(II+4, KOD(I,5)) =  Q(3)-Q(1);

%   GG(II+1, KOD(I,6)) = 1-RS(I,3)       ;
%   GG(II+3, KOD(I,6)) = RS(I,3)         ;
%   GG(II+4, KOD(I,6)) = RS(I,3) -RS(I,1);

%   GG(II+2, KOD(I,7)) =  PT(I,2)        ;
%   GG(II+4, KOD(I,7)) = -PT(I,1)-PT(I,2);
%   GG(II+5, KOD(I,7)) = -PT(I,2)        ;

%   GG(II+2, KOD(I,8)) = -Q(2)           ;
%   GG(II+4, KOD(I,8)) =  Q(2)-Q(1)      ;
%   GG(II+5, KOD(I,8)) =  Q(2)           ;

%   GG(II+2, KOD(I,9)) = 1-RS(I,2)       ;
%   GG(II+4, KOD(I,9)) = RS(I,2)-RS(I,1) ;
%   GG(II+5, KOD(I,9)) = RS(I,2)         ;
% end

% XYZ = [  B(2)*GG(1,:)+B(3)*GG(2,:)   ;
%        2*B(2)*GG(3,:)+B(3)*GG(4,:)   ;
%          B(2)*GG(4,:)+2*B(3)*GG(5,:) ;
%         -C(2)*GG(6,:)-C(3)*GG(7,:);
% QQ(5,:) = -2*C(2)*GG(8,:) - C(3)*GG(9,:) ;
% QQ(6,:) = -C(2)*GG(9,:) -2*C(3)*GG(10,:) ;

% QQ(7,:) =  C(2)*GG(1,:) + C(3)*GG(2,:) ...
%           -B(2)*GG(6,:) - B(3)*GG(7,:) ;

% QQ(8,:) = 2*C(2)*GG(3,:) + C(3)*GG(4,:) ...
%          -2*B(2)*GG(8,:) - B(3)*GG(9,:) ;

% QQ(9,:) = C(2)*GG(4,:) + 2*C(3)*GG(5,:) ...
%          -B(2)*GG(9,:) - 2*B(3)*GG(10,:) ;

% K = QQ' * DD * QQ ;



% implementation of the DKT plate triangle element based on https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620180711

x12 = x1 - x2 ;   x23 = x2 - x3 ;   x31 = x3 - x1 ;
y12 = y1 - y2 ;   y23 = y2 - y3 ;   y31 = y3 - y1 ;

mat_cross = [ x23  x31; y23  y31; 0  0 ] ;
Area = norm( cross( mat_cross(:,1), mat_cross(:,2)) )*.5 ;

l212 = x12^2 + y12^2 ;
l223 = x23^2 + y23^2 ;
l231 = x31^2 + y31^2 ;

p4 = -6*x23/l223 ;
p5 = -6*x3/l231  ;
p6 = -6*x12/l212 ;

t4 = -6*y23/l223 ;
t5 = -6*y3/l231  ;

q4 = 3*x23*y23/l223 ;
q5 = 3*x3*y3/l231 ;

r4 = 3*y23^2/l223 ;
r5 = 3*y31^2/l231 ;


alphaMat = [ ...
  y3*p6        0          -4*y3       -y3*p6         0          -2*y3       0            0                0             ;
 -y3*p6        0           2*y3        y3*p6         0          4*y3        0            0                0             ;
  y3*p5       -y3*q5       y3*(2-r5)   y3*p4         y3*q4      y3*(r4-2)  -y3*(p4+p5)   y3*(q4-q5)       y3*(r4-r5)    ;
 -x2*t5        x23+x2*r5  -x2*q5       0             x3         0           x2*t5         x2*(r5-1)       -x2*q5        ;
  0            x23         0           x2*t4         x3+x2*r4  -x2*q4      -x2*t4         x2*(r4-1)       -x2*q4        ;
  x23*t5       x23*(1-r5)  x23*q5     -x3*t4         x3*(1-r4)  x3*q4      -x23*t5+x3*t4 -x23*r5-x3*r4-x2  x3*q4+x23*q5 ;
 -x3*p6-x2*p5  x2*q5+y3   -4*x23+x2*r5 x3*p6        -y3         2*x3        x2*p5         x2*q5            (r5-2)*x2    ;
 -x23*p6       y3          2*x23       x23*p6+x2*p4 -y3+x2*q4  -4*x3+x2*r4 -x2*p4         x2*q4            (r4-2)*x2    ;
 %
  x23*p5+y3*t5  -x23*q5+(1-r5)*y3   (2-r5)*x23+y3*q5   -x3*p4+y3*t4  (r4-1)*y3-x3*q4 (2-r4)*x3-y3*q4 ...
  %
  -x23*p5+x3*p4-(t4+t5)*y3    -x23*q5-x3*q4+(r4-r5)*y3  -x23*r5-x3*r4+4*x2+(q5-q4)*y3 ] ;

% isotropic
Db = thickness^3/12 * E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)*.5 ] ;

R = ones(3,3)+eye(3) ;

DL =  zeros(9,9);
i = 1; j = 1;  DL(nodes2dofs(i,3),nodes2dofs(j,3)) = 1/24* Db(i,j)*R ;
i = 1; j = 2;  DL(nodes2dofs(i,3),nodes2dofs(j,3)) = 1/24* Db(i,j)*R ;
i = 2; j = 1;  DL(nodes2dofs(i,3),nodes2dofs(j,3)) = 1/24* Db(i,j)*R ;
i = 2; j = 2;  DL(nodes2dofs(i,3),nodes2dofs(j,3)) = 1/24* Db(i,j)*R ;
i = 3; j = 3;  DL(nodes2dofs(i,3),nodes2dofs(j,3)) = 1/24* Db(i,j)*R ;

K = 1/(2*Area) * alphaMat' * DL * alphaMat ;

fint = K * elemDisps ;

fs = {fint};  ks = {K};

% fout = zeros(18,1);
% fout(dofs_sort) = fint ;

% Kout = zeros(18,18) ;
% Kout(dofs_sort,dofs_sort) = K ;

% fs = {fout};  ks = {Kout};