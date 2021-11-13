% ======================================================================
% Truss-beam Joint example
% ======================================================================

close all, clear all

dirOnsas = [ pwd '/../../src' ] ; % set ONSAS.m directory
addpath( dirOnsas ); % add ONSAS directory to path
problemName = 'TurssBeamJoint_NR' ; 

% ----------------------------------------------------------------------
%Geometric properties truss:
Et = 1e9 ; dt = .05; At = pi*dt^2/4 ;  Lt = 1 ; nut = 0.3 ;  rhot = 1 ; 
%Geometric properties beam
Eb = Et/3 ;db = 5*dt ; Ab = pi*db^2/4 ;  Lb = .5 ; nub = 0.3 ;  rhob = 1 ; 
Ib = pi*db^4/64 ;
% ----------------------------------------------------------------------
% MELCS parameters

% Materials
materialsParams = {[ rhot 1 Et nut ];[rhob 1 Eb nub ]} ;

% Elements
elementsParams  = {1;2;3} ;

% Loads
loadsParams     = { [ 1 1   0 0 0 0 +1 0] } ;

% Cross-Sections
crossSecsParams = {[3 dt] ; [3 db]} ;	

% Springs
springsParams   = { [ inf  inf    inf  inf  inf   inf  ] ; ...
                    [  0   inf    inf   0    0     inf ] ; ...
                    [ inf  inf    inf   inf    inf   0 ] } ;

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


Nodes = [ 	(0:(NelemB))'*Lb/NelemB  zeros(NelemB+1,1)          zeros(NelemB+1,1) 	   ;	 
			Lb*ones(NelemT,1) 		 zeros(NelemT,1)         -(1:(NelemT))'*Lt/NelemT ] ;
													   
				    %Material             %Element                   %Loads                 %CrossSection             %Springs 					 
auxConecElem  = [ [(ones(NelemB,1)*2 )    (ones(NelemB,1)*3)      (zeros(NelemB,1))         (ones(NelemB,1)*2)      (zeros(NelemB,1)) ...  %ElemNodes...
                 (1:(NelemB))'                 (2:NelemB+1)'                  zeros(NelemB,2)];
                  [(ones(NelemT,1)*1 )    (ones(NelemT,1)*2)      (zeros(NelemT,1))         (ones(NelemT,1)*1)      (zeros(NelemT,1)) ...  %ElemNodes...
                 (NelemB+1:NelemB+NelemT)'	 	(NelemB+2:NelemB+NelemT+1)'   zeros(NelemT,2)] ] ;  
                        %Material             %Element               ma    %Loads                 %CrossSection             %Springs 			%Nodes  
auxConecNodes = [           0                   1                           0                           0                   1                   1   
                            0                   1                           1                           0                   2                NnodesB 
                            0                   1                           0                           0                   3                 Nnodes ];                                                                                                                                 

stop
%Build Conec Cell
Conec = cell (Nelem+size(auxConecNodes,1),1);
for i = 1:size(auxConecNodes,1)
    Conec{i,1}=auxConecNodes(i,:) ;
end
for i =  1:Nelem
    Conec{size(auxConecNodes,1)+i,1} =auxConecElem(i,:);
end

controlDofs      = [ NnodesB 5 1 ] ; % [ node nodaldof scalefactor ]
plotParamsVector = [ 3 ] ;
reportBoolean    = 1     ;
storeBoolean     = 1     ;

% analysis parameters
stopTolDeltau    = 1.0e-12 ;    stopTolForces    = 1.0e-12;
targetLoadFactr  = 1e4 ;    nLoadSteps       = 10      ;
stopTolIts       = 10     ;
global flagOutputMatrices;
flagOutputMatrices= 1 ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 

%Analytic Solution
analyticSolFlag = 2 
analyticCheckTolerance = 1e-5 ;
analyticFunc   = @(w)(Et*At/Lt+3*Eb*Ib/Lb^3)*w;
beamTruss_Ratio= Et*At/Lt/(3*Eb*Ib/Lb^3)

% run ONSAS                        
ONSAS


controlDispsNR = controlDisps ;
loadFactorsNR  = loadFactors ;

% ----------------------------------------------------------------------

problemName = 'TurssBeamJoint_AL' ; 

%Build Conec Cell
Conec = cell (Nelem+size(auxConecNodes,1),1);
for i = 1:size(auxConecNodes,1)
    Conec{i,1}=auxConecNodes(i,:) ;
end
for i =  1:Nelem
    Conec{size(auxConecNodes,1)+i,1} =auxConecElem(i,:);
end

% arc length params
targetLoadFactrNRAL   = targetLoadFactr  ;
incremArcLen          = controlDispsNR(end)- controlDispsNR(end-1) ;%Same DeltaU step
numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactrNRAL nLoadSteps incremArcLen ] ; 
                        
ONSAS
controlDispsAL = controlDisps ;
loadFactorsAL  = loadFactors ;



% ----------------------------------------------------------------------
% --- plots ---
lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;

figure
plot( controlDispsNR, analyticFunc(controlDispsNR),'b-x' , 'linewidth', lw,'markersize',ms)
hold on, grid on
plot( controlDispsAL, loadFactorsAL,'r-s' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNR, loadFactorsNR,'k-o' , 'linewidth', lw,'markersize',ms )

labx = xlabel('Displacement');   laby = ylabel('\lambda') ;
% legend('analytic','NRAL-DXF','NR','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;



