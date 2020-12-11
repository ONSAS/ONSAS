% ======================================================================
% Truss-beam Joint example

close all, clear all %#ok

dirOnsas = [ pwd '/../..' ] ; % set ONSAS.m directory
addpath( dirOnsas ); % add ONSAS directory to path

problemName = 'cantileverPendulum' ; 
% ----------------------------------------------------------------------
%% Geometric Properties
%Geometric properties truss:
Et = 10e11 ; At = .1 ; dt = sqrt(4*At/pi);   Lt = 3.0443 ; nut = 0.3 ;  rhot = 65.6965 ; 
%Geometric properties beam
Eb = Et/300000*7 ;db = dt ; Ab = pi*db^2/4 ;  Lb = Lt ; nub = 0.3 ;  rhob = rhot ; 
Ib = pi*db^4/64 ;
% ----------------------------------------------------------------------
%% MELCS parameters

% Materials
materialsParams = {[ rhot 3 Et nut ]; [rhob 3 Eb nub ]} ;

% Elements
elementsParams  = {1; [2 1]; 3} ;

% Loads
loadsParams     = { [ 1 1   0 0 0 0 1 0] } ;

% Cross-Sections
crossSecsParams = {[3 dt] ; [3 db]} ;	

% Springs
springsParams   = { [ inf  inf    inf  inf  inf   inf  ]} ;

%Nodes and Conec matrix:
%-----------------------

%Truss elements
NelemT  = 1;
NnodesT = NelemT+1;

%Beam elements
NelemB  = 10;
NnodesB = NelemB+1;

%Total elements
Nelem = NelemT+NelemB	;
Nnodes = Nelem+1;

Nodes = [ 	(0:(NelemB))'*Lb/NelemB        zeros(NelemB+1,1)     zeros(NelemB+1,1) ;	 
		(Lb+(1:(NelemT))'*Lt/NelemT)  zeros(NelemT,1)       zeros(NelemT,1)] ;
													   
		    %Material             Element                   Loads                 CrossSection     Springs 					 
auxConecElem  = [ [(ones(NelemB,1)*2 )    (ones(NelemB,1)*3)      (zeros(NelemB,1))         (ones(NelemB,1)*2)   (zeros(NelemB,1)) ...  %ElemNodes...
                 (1:(NelemB))'                 (2:NelemB+1)'                  zeros(NelemB,2)];
                  [(ones(NelemT,1)*1 )    (ones(NelemT,1)*2)      (zeros(NelemT,1))         (ones(NelemT,1)*1)   (zeros(NelemT,1)) ...  %ElemNodes...
                 (NelemB+1:NelemB+NelemT)'	 	(NelemB+2:NelemB+NelemT+1)'   zeros(NelemT,2)] ] ;  
                        %Material             %Element               ma    %Loads                 %CrossSection             %Springs 			%Nodes  
auxConecNodes = [           0                   1                           0                           0                   1                   1   
                            0                   1                           1                           0                   0                 Nnodes ];

% Conec Cell
Conec = cell (Nelem+size(auxConecNodes,1),1);
for i = 1:size(auxConecNodes,1)
    Conec{i,1}=auxConecNodes(i,:) ;
end
for i =  1:Nelem
    Conec{size(auxConecNodes,1)+i,1} =auxConecElem(i,:);
end

% ----------------------------------------------------------------------
%% Loads parameters
booleanSelfWeightZ = 1 ;
loadsParams   = {[ 1 1   0 0 0 0 -1 0 ]} ;
massPendulum  = 5 ;
gravity       = 9.8;
loadFactorsFunc = @(t)0* massPendulum*gravity;
% ----------------------------------------------------------------------
%% analysis parameters
% method
timeIncr   =  0.05    ;
finalTime  =  4.13*8 ;
alphaHHT   = -0.05 ;

% tolerances
stopTolDeltau = 1e-5 ; 
stopTolForces = 1e-5 ;
stopTolIts    = 30   ;
% ------------------------------------
% analysis parameters
numericalMethodParams = [ 4 ...
  timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;

%Output params
controlDofs      = [ NnodesT 5 -1 ] ; % [ node nodaldof scalefactor ]
plotParamsVector = [ 3 ] ;
reportBoolean    = 1     ;
storeBoolean     = 1     ;
% ----------------------------------------------------------------------

% Run ONSAS
ONSAS

% ----------------------------------------------------------------------

%% --- plots ---
lw = 2.0 ; ms = 15 ; plotfontsize = 15 ;

figure
time = 0:timeIncr:finalTime+timeIncr;
plot( time,controlDisps,'b-' , 'linewidth', lw,'markersize',ms)


labx = xlabel('time (s)');   laby = ylabel('Displacement (m)') ;
% legend('analytic','NRAL-DXF','NR','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
