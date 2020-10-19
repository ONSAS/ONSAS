% ======================================================================
% Truss-beam Joint example

close all, clear all %#ok
problemName = 'cantileverPendulum' ; 
% ----------------------------------------------------------------------
%% Geometric Properties
%Geometric properties truss:
Et = 1e9 ; dt = .05; At = pi*dt^2/4 ;  Lt = 1 ; nut = 0.3 ;  rhot = 1 ; 
%Geometric properties beam
Eb = Et/3 ;db = 5*dt ; Ab = pi*db^2/4 ;  Lb = .5 ; nub = 0.3 ;  rhob = 1 ; 
Ib = pi*db^4/64 ;
% ----------------------------------------------------------------------
%% MELCS parameters

% Materials
materialsParams = {[ rhot 1 Et nut ];[rhob 1 Eb nub ]} ;

% Elements
elementsParams  = {1;[2 1];3} ;

% Loads
loadsParams     = { [ 1 1   0 0 0 0 +1 0] } ;

% Cross-Sections
crossSecsParams = {[3 dt] ; [3 db]} ;	

% Springs
springsParams   = { [ inf  inf    inf  inf  inf   inf  ]} ;

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


Nodes = [ 	(0:(NelemB))'*Lb/NelemB        zeros(NelemB+1,1)     zeros(NelemB+1,1) ;	 
			(Lb+(1:(NelemT))'*Lt/NelemT)  zeros(NelemT,1)       zeros(NelemT,1)] ;
													   
				    %Material             %Element                   %Loads                 %CrossSection             %Springs 					 
auxConecElem  = [ [(ones(NelemB,1)*2 )    (ones(NelemB,1)*3)      (zeros(NelemB,1))         (ones(NelemB,1)*2)      (zeros(NelemB,1)) ...  %ElemNodes...
                 (1:(NelemB))'                 (2:NelemB+1)'                  zeros(NelemB,2)];
                  [(ones(NelemT,1)*1 )    (ones(NelemT,1)*2)      (zeros(NelemT,1))         (ones(NelemT,1)*1)      (zeros(NelemT,1)) ...  %ElemNodes...
                 (NelemB+1:NelemB+NelemT)'	 	(NelemB+2:NelemB+NelemT+1)'   zeros(NelemT,2)] ] ;  
                        %Material             %Element               ma    %Loads                 %CrossSection             %Springs 			%Nodes  
auxConecNodes = [           0                   1                           0                           0                   1                   1   
                            0                   1                           1                           0                   0                 Nnodes ];                                                                                                                                 


%Build Conec Cell
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
loadsParams   = {[ 1 1   0 0 0 0 1 0 ]} ;
massPendulum  = 10 ;
gravity       = 9.8;
loadFactorsFunc = @(t) -massPendulum*gravity;
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
controlDofs      = [ NnodesB 5 1 ] ; % [ node nodaldof scalefactor ]
plotParamsVector = [ 3 ] ;
reportBoolean    = 1     ;
storeBoolean     = 1     ;
% ----------------------------------------------------------------------
%% RunONSAS
beamTruss_Ratio= Et*At/Lt/(3*Eb*Ib/Lb^3)
                        
run( [ pwd '/../ONSAS.m' ] ) ;

controlDispsNR = controlDisps ;
loadFactorsNR  = loadFactors ;


% ----------------------------------------------------------------------
%% --- plots ---
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
