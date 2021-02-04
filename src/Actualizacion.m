%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%				Actualización de variables	 	 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Dtp1kp1,Ddottp1kp1,Ddotdottp1kp1] = Actualizacion (Dt,Dtp1k,Ddott,Ddottp1k,Ddotdott,Ddotdottp1k,DeltaDtp1,NewmarkParams,numericalMethodParams,Nnodes,IterDisp,IterTime)


%Newmark paramas (estan al reves que en ANLE) (beta =1/2,Alpha=1/4)
beta   = NewmarkParams(1);
alpha  = NewmarkParams(2);
DeltaT = numericalMethodParams(4);

%---------------- Initialize Disp, Ang Vel and Acel iter k  -------------------

Dtp1kp1		  	 = zeros(6*Nnodes,1) ;
Ddottp1kp1	  	 = zeros(6*Nnodes,1) ;
Ddotdottp1kp1    = zeros(6*Nnodes,1) ;
DeltaThetatp1kp1 = zeros(3*Nnodes,1);

%---------------- Separacion Ang-Disp  -------------------
DeltaUtp1k 	                = DeltaDtp1(1:2:end) ;
DeltaThetatp1k 				= DeltaDtp1(2:2:end) ; % esta al reves la convencion ? lo cambie por el inpou NOdal MOment


%Displacements t
Ut		 = Dt (1:2:end)		 ;
Udott	 = Ddott (1:2:end)	 ;
Udotdott = Ddotdott (1:2:end);

%Angles t
Thetat   = Dt (2:2:end)		 ;
Wdott	 = Ddott (2:2:end)	 ;
Wdotdott = Ddotdott (2:2:end);

%Displacements t+Delta_t paso k

%Displacements 
Utp1k		  = Dtp1k       (1:2:end) ;
Udottp1k	  = Ddottp1k    (1:2:end) ;
Udotdottp1k   = Ddotdottp1k (1:2:end) ;

%Angles 
Thetatp1k   = Dtp1k       (2:2:end)	;
Wdottp1k    = Ddottp1k    (2:2:end)	;
Wdotdottp1k = Ddotdottp1k (2:2:end) ;



%---------------- Actualizacion Disp  -------------------
Utp1kp1		   = Utp1k + DeltaUtp1k;

Udotdottp1kp1  = 1/(beta*(DeltaT^2))*(Utp1kp1-Ut) - 1/(beta*DeltaT)*Udott...
						- (1/(2*beta))*(1-2*beta) * Udotdott		;% Ec 2.10

Udotp1kp1 	   = alpha/(beta*(DeltaT))*(Utp1kp1-Ut) + (1-alpha/beta) * Udott...
						+ (1 - alpha/(2*beta))* Udotdott		    ; %Ec 2.11

%---------------- Actualizacion Angulo  -------------------

Thetatp1kp1 = zeros(size (DeltaThetatp1k,1),1);

%La funcion skew es para tres entonces hace nodo a nodo:

for Node=1:Nnodes
	
	index1 = 3*Node-2;
	index2 = 3*Node  ;

	DeltaThetatp1i = DeltaThetatp1k(index1:index2)	; % elige los incrementos angulos correspondiente al alemento
	Thetatp1ki     = Thetatp1k (index1:index2)      ; % elige el valor de angulo actual de theta
	

	%~ Thetatp1ki
	Lamdai 		   = expon(DeltaThetatp1i) * expon(Thetatp1ki)  ; %calcula lamdai EC 99

	%Arreglo del orden de la base [v(1),-v(3),v(2)]
	% Calculo de los angulos
	Thetatp1kp1i                    = logar(Lamdai);%Inverso Ec 96
	%~ Thetatp1kp1i 					= [Thetatp1kp1i(1);-Thetatp1kp1i(3);Thetatp1kp1i(2)] %Convencion:
	Thetatp1kp1  (index1:index2)    = Thetatp1kp1i;%Coloco en el vector total
	
	Lamdai                          = expon(Thetatp1kp1i);

	% Calculo de las velocidades
	Wdottp1kp1i						=  Lamdai * [alpha/(beta*DeltaT)*Thetatp1kp1i + (beta-alpha)/beta*Wdott(index1: index2) + (beta - 0.5*alpha)/beta*Wdotdott(index1: index2) ];%Ec 100
	%~ Wdottp1kp1i 					=  [-Wdottp1kp1i(1),-Wdottp1kp1i(2),-Wdottp1kp1i(3)];
	Wdottp1kp1   (index1:index2)    =  Wdottp1kp1i ;
	
	% Calculo de las acleraciones
	Wdotdottp1kp1i 					= Lamdai * [1/(beta*DeltaT^2)*Thetatp1kp1i - 1/(beta*DeltaT)*Wdott(index1: index2) - (.5 - beta)/beta*Wdotdott(index1: index2) ] ; %Ec 101
	%~ Wdottp1kp1i 					= [-Wdottp1kp1i (1),-Wdottp1kp1i (2),-Wdottp1kp1i (3)];
	Wdotdottp1kp1(index1:index2) 	= Wdotdottp1kp1i;
end	
	
%~ Utp1kp1
%~ Udotp1kp1
%~ Udotdottp1kp1
%~ fprintf('Angulares......:\n')
%~ Thetatp1kp1
%~ Wdottp1kp1
%~ Wdotdottp1kp1


%------ Comptactar desplazamientos velcodiadades y aceleraciones ----------
Dtp1kp1(1:2:end) = Utp1kp1;
Dtp1kp1(2:2:end) = Thetatp1kp1;

Ddottp1kp1(1:2:end) = Udotp1kp1;
Ddottp1kp1(2:2:end) = Wdottp1kp1;

Ddotdottp1kp1(1:2:end) = Udotdottp1kp1;
Ddotdottp1kp1(2:2:end) = Wdotdottp1kp1;
end


