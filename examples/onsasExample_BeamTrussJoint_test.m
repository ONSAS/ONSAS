% ======================================================================
% Truss-beam Joint example

close all, clear all %#ok
problemName = 'TurssBeamJoint' ; 
% ----------------------------------------------------------------------
%Geometric properties truss:
Et = 210e5 ; dt = 0.03; At = pi*dt^2/4 ;  Lt = 1 ; nut = 0.3 ;  rhot = 8000 ; 
%Geometric properties beam
Eb = 210e9 ; db = 0.1; Ab = pi*db^2/4 ;  Lb = 10 ; nub = 0.3 ;  rhob = 1000 ; 

% ----------------------------------------------------------------------
% MELCS parameters

% Materials
materialsParams = {[ rhot 1 Et nut ];[rhob 1 Eb nub ]} ;

% Elements
elementsParams  = {1;2;3} ;

% Loads
loadsParams     = { [ 1 1   0 0 0 0 -1 0] } ;

% Cross-Sections
crossSecsParams = {[3 dt] ; [3 db]} ;	

% Springs
springsParams   = { [ inf  inf  inf  inf  inf   inf ] ; ...
                    [ inf    0  inf  0    inf   0 ] } ;

%Nodes and Conec matrix:
%-----------------------

%Truss elements
NelemT  = 1;
NnodesT = NelemT+1;

%Bem elements
NelemB  = 1;
NnodesB = NelemB+1;

%Total elements
Nelem = NelemT+NelemB	;
Nnodes = Nelem+1;


Nodes = [ 	(0:(NelemT))'*Lb/NelemT  zeros(NnodesT,1)          zeros(NnodesT,1) 	   ;	 
			Lb*ones(NelemB,1) 		 zeros(NelemB,1)         -(1:(NelemB))'*Lt/NelemB ] ;
													   
				    %Material             %Element                   %Loads                 %CrossSection             %Springs 					 
auxConecElem  = [ [(ones(NelemB,1)*2 )    (ones(NelemB,1)*3)      (ones(NelemB,1))         (ones(NelemB,1)*2)      (zeros(NelemB,1)) ...  %ElemNodes...
                 (1:(NelemB))'                 (2:NelemB+1)'                  zeros(NelemB,2)];
                  [(ones(NelemT,1)*1 )    (ones(NelemT,1)*2)      (ones(NelemT,1))         (ones(NelemT,1)*1)      (zeros(NelemT,1)) ...  %ElemNodes...
                 (NelemB+1:NelemB+NelemT)'	 	(NelemB+2:NelemB+NelemT+1)'   zeros(NelemB,2)] ] ;  
                        %Material             %Element                   %Loads                 %CrossSection             %Springs 			%Nodes  
auxConecNodes = [           0                   1                           0                           0                   1                   1   
                            0                   1                           1                           0                   0                NnodesB 
                            0                   1                           0                           0                   2                 Nnodes ];                                                                                                                                 


%Build Conec Cell
Conec = cell (Nelem+Nnodes,1);
for i = 1:Nnodes
    Conec{i,1}=auxConecNodes(i,:) ;
end
for i =  1:Nelem
    Conec{Nnodes+i,1} =auxConecElem(i,:);
end

controlDofs      = [ 2 5 -1 ] ; % [ node nodaldof scalefactor ]
plotParamsVector = [ 3 ] ;
reportBoolean    = 1     ;


% analysis parameters
stopTolDeltau    = 1.0e-8 ;    stopTolForces    = 1.0e-8 ;
targetLoadFactr  = 2.0e10  ;    nLoadSteps       = 20      ;
stopTolIts       = 30     ;


numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 

run( [ pwd '/../ONSAS.m' ] ) ;

