clear all
close all
clc
%
nelems = 4 ;
states = 3 ;

%
state = cell(nelems,states) ;
%
state(:,1) = zeros(1,3)
state(:,2) = zeros(1,3)
state(:,3) = zeros(1,3)



%%%
elem = 2 ;

stress = state(:,1)
strain = state(:,1)
acum_p_strain = state(:,3) 

stateVector = [ stress{elem} ; strain{elem} ; acum_p_strain{elem} ]

stateVector(1)

% ======================================================================
E =  1000000
nu =  0.30000
G=E/(2*(1+nu))
K = E/(3*(1-2*nu)) +G/3; % Bulk modulus

%~ strain = [   0.0103075; 0.0031424; -0.0385666]

De = E * (1-nu) / ( (1+nu)*(1-2*nu) ) * ...
											 [ 1          nu/(1-nu)  0                   ; ...
												 nu/(1-nu)  1          0                   ; ...
												 0          0          (1-2*nu)/(2*(1-nu)) ] ;

%~ strain= [47.375   39.933  -10.014]
						 

%%
strain=[
  -0.0105703
   0.0216760
  -0.0022620
]

Stress = De*strain ;

StressTensor = [Stress(1) Stress(3) ; Stress(3) Stress(2)]

StrainTensor = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)]

HydStressTensor = 1/2*trace(StressTensor)*eye(2) 
DevStressTensor = StressTensor - HydStressTensor 

HydStrainTensor = 1/2*trace(StrainTensor)*eye(2) 
DevStrainTensor = StrainTensor - HydStrainTensor 

HydStrainVec = [ HydStrainTensor(1,1) HydStrainTensor(2,2) 0 ]
DevStrainVec = [ DevStrainTensor(1,1) DevStrainTensor(2,2) 2*DevStrainTensor(1,2) ]

st = De*(HydStrainVec+DevStrainVec)'

hyd_st = De*HydStrainVec'
dev_st = De*DevStrainVec'

hyd_st = K*sum(HydStrainVec)*eye(2)

DevStrainVec = [ DevStrainTensor(1,1) DevStrainTensor(2,2) DevStrainTensor(1,2) ]
dev_st = 2*G*DevStrainVec


%~ a = [ 1 2 0 ; 2 5 0 ; 0 0 (1+5)/2 ]

%~ hyd = trace(a)/3*eye(3)
%~ dev = a - hyd

%~ aa = [ 1 2 ; 2 5 ]
%~ hyd = trace(aa)/2*eye(2)
%~ dev = aa - hyd

Is = eye(3) ;
Is(3,3) = 1/2 ;

Id = Is - 1/3*Is*trace(Is) ;

De2 = 2*G*Id + K*eye(3) ;
