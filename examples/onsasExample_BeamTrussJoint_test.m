% ======================================================================
% Truss-beam Joint example

close all, clear all %#ok
problemName = 'TurssBeamJoint_NR' ; 
% ----------------------------------------------------------------------
%Geometric properties truss:
Et = 210e5 ; dt = 0.03; At = pi*dt^2/4 ;  Lt = 1 ; nut = 0.3 ;  rhot = 8000 ; 
%Geometric properties beam
Eb = 210e9 ; db = 0.1; Ab = pi*db^2/4 ;  Lb = 10 ; nub = 0.3 ;  rhob = 1000 ; 
Ib = pi*db^4/64 ;
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
                    [ inf  inf  inf  inf   0     0 ] ; ...
                    [ inf    0  inf  0    inf   0 ] } ;

%Nodes and Conec matrix:
%-----------------------

%Truss elements
NelemT  = 1;
NnodesT = NelemT+1;

%Bem elements
NelemB  = 2;
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


%Build Conec Cell
Conec = cell (Nelem+size(auxConecNodes,1),1);
for i = 1:size(auxConecNodes,1)
    Conec{i,1}=auxConecNodes(i,:) ;
end
for i =  1:Nelem
    Conec{size(auxConecNodes,1)+i,1} =auxConecElem(i,:);
end

controlDofs      = [ NnodesB 5 -1 ] ; % [ node nodaldof scalefactor ]
plotParamsVector = [ 3 ] ;
reportBoolean    = 1     ;
storeBoolean     = 1     ;

% analysis parameters
stopTolDeltau    = 1.0e-8 ;    stopTolForces    = 1.0e-8 ;
targetLoadFactr  = 2.0e11  ;    nLoadSteps       = 20      ;
stopTolIts       = 30     ;


numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 

run( [ pwd '/../ONSAS.m' ] ) ;


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
targetLoadFactrNRAL   = 2e11  ;
incremArcLen          = 1e-2  ;
numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactrNRAL nLoadSteps incremArcLen ] ; 
                        
run( [ pwd '/../ONSAS.m' ] ) ;
controlDispsNRAL = controlDisps ;
loadFactorsNRAL  = loadFactors ;

analyticFunc =@(w) (Et*At/Lt+3*Eb*Ib/Lb^3)*w;

% ----------------------------------------------------------------------
% --- plots ---
lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;

figure
plot( controlDispsNRAL, analyticFunc(controlDispsNRAL),'b-x' , 'linewidth', lw,'markersize',ms)
hold on, grid on
% plot( controlDispsNRAL, loadFactorsNRAL,'r-s' , 'linewidth', lw,'markersize',ms )
% plot( controlDispsNR, loadFactorsNR,'k-o' , 'linewidth', lw,'markersize',ms )

labx = xlabel('Displacement');   laby = ylabel('\lambda') ;
% legend('analytic','NRAL-DXF','NR','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;



