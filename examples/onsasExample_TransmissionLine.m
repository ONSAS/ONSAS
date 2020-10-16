% ------------------------------------------------------------------------------
% example TransmissionLine
%
% ------------------------------------------------------------------------------

clear all, close all
inputONSASversion = '0.1.10';
%acdir = pwd ; 
problemName = 'TransmissionLine_Selfweight' ;



%Geometric properties cable:
%-------------------- 
dc  = 0.1 ; lc  = 406.5;
% Ac  = pi*dc^2/4 ;
% Iyc = pi*dc^4/64;
% Izc = Iyc        ; 
% Jc  = 2*Iyc     ;

%Geometric properties isolator chain:
%-------------------- 
da  = 0.1 ; la  = 5;
% Aa  = pi*da^2/4  ;
% Iya = pi*da^4/64 ;
% Iza = Iya        ; 
% Ja  = 2*Iya      ;

%Material properties cable :
%--------------------
Ec   = 70e9    ;
nuc  = 0.3      ; 
rhoc = 2700	    ;
Gc	= Ec/(2*(1+nuc));	

%Material properties isolator chain :
%--------------------
Ea   = 72e9    ;
nua  = 0.3      ; 
rhoa = 2555	    ;
Ga	= Ea/(2*(1+nua));	

%Elem CrosSec and elem Params Cells:
%--------------------
%1 cable 2 isolator chain
materialsParams = {[rhoc 1 Ec nuc ] ;[rhoa 1 Ea nua ]} ;

crossSecsParams = {[2 dc dc] ; [2 da da]} ;	
% crossSecsParams = {[3 dc] ; [3 da]} ;	
elementsParams  = { 1; 3} ;

%Nodes and Conec matrix:
%-----------------------

%Nelem Cable
NelemC  = 20;
NnodesC = NelemC+1;

%Nelem Aisladora
NelemA  = 2;
NnodesA = NelemA+1;

%Nelem Tottales
Nelem = NelemC+NelemA*2	;
Nnodes = Nelem+1;


Nodes = [ 	  zeros(NelemA,1) 		 zeros(NelemA,1)        ((NelemA):-1:1)'*la/NelemA ;
			(0:(NelemC))'*lc/NelemC  zeros(NnodesC,1)          zeros(NnodesC,1) 	   ;	 
			lc*ones(NelemA,1) 		 zeros(NelemA,1)         (1:(NelemA))'*la/NelemA ] ;
													   
				   %Material             %Element                   %Loads                 %CrossSection             %Springs 					 
auxConecElem  =[ [(ones(NelemA,1)*2 )    (ones(NelemA,1)*2)      (ones(NelemA,1))         (ones(NelemA,1)*2)      (zeros(NelemA,1)) ...  %ElemNodes...
                (1:(NelemA))'                               (2:(NnodesA))'                              zeros(NelemA,2)];  
                 [(ones(NelemC,1)*1 )    (ones(NelemC,1)*2)      (ones(NelemC,1))         (ones(NelemC,1)*1)      (zeros(NelemC,1)) ...  %ElemNodes...
                (NelemA+1:(NelemA+NelemC))'                 (NelemA+2:NelemA+1+NelemC)'                  zeros(NelemC,2)];
                 [(ones(NelemA,1)*2 )    (ones(NelemA,1)*2)      (ones(NelemA,1))         (ones(NelemA,1)*2)      (zeros(NelemA,1)) ...  %ElemNodes...
                (NelemA+NelemC+1:NelemA+NelemC+NelemA)'	 	(NelemA+NelemC+2:NelemA+NelemC+NelemA+1)']   zeros(NelemA,2)  ] ;  
  
auxConecNodes =[ (zeros(Nnodes,1)*1 )    (ones(Nnodes,1))         (ones(Nnodes,1))         (zeros(Nnodes,1))        (zeros(Nnodes,1))    (1:Nnodes)'] ;

% Fixed nodes
auxConecNodes(1,5)              = 1 ;
auxConecNodes(end,5)            = 1 ;
auxConecNodes(NnodesA,5)        = 2 ;
auxConecNodes(Nnodes-NnodesA+1,5) = 2 ;


%Build Conec Cell
Conec = cell (Nelem+Nnodes,1);
for i = 1:Nnodes
    Conec{i,1}=auxConecNodes(i,:) ;
end
for i =  1:Nelem
    Conec{Nnodes+i,1} =auxConecElem(i,:);
end
 
%Suports---------------------
springsParams    = {[ inf  0  inf  0  inf  0 ];[ inf  0  0  inf  0  0 ]} ;

%Parameters and tolerances:
%-----------------------
% times
finalTime  = 200    ;
timeIncr   =  finalTime/400    ;

nLoadSteps = finalTime/timeIncr ;

% tolerances
stopTolDeltau = 1e-3        ; 
stopTolForces = 1e-7        ;
stopTolIts    = 30          ;

%method params
%Newmark
%~ DeltaNW    =  0.5               ;
%~ AlphaNW    =  0.25              ;
%~ numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

%HHT
alphaHHT = 0 ;
alphaHHT = -0.05 ;
numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;
    
%Loads and supports:
%Loads---------------------
booleanSelfWeightZ = 1 ;
loadsParams   = {[ 1 1   0 0 1 0 0 0 ]} ;

loadFactorsFunc = @(t)0*t;
% Damping
nodalDispDamping = 1;
%Booleans control and plot:
%-----------------------
%controlDof
controlDofs = [ Nelem/2 5 1 ] ;
%Plots
plotParamsVector = [ 3 150 ];
printFlag = 0 ;
%Booleans
storeBoolean = 1  ;
reportBoolean = 0 ;

%RunONSAS
run( [ pwd '/../ONSAS.m' ] ) ;