## Modelo elastico
## Author: NGT

##Rev0: 9/2/19 Unifico módulos. Uso de plantilla aero3D_ala_rev4
##Rev1: 9/2/19 Agrego subiteraciones en estructura para trabajar con funciones de forma de alta frecuencia
##Rev2: ./2/19 
##Rev3: 20/2/19 Uso integrador Newmark-beta para estructura en vez de subiteraciones (m�s r�pido y estable)

warning ("off", "Octave:broadcast");
warning ("off", "Octave:possible-matlab-short-circuit-operator");

clear all
close all
clc

restart=0;

if restart==0
###########################################################
#Parametros para integradores temporales
deltaT=1; #Siempre 1
#Parametros simulacion 
rho=1.225; #densidad del fluido, aire SI
wy=2 #velocidad rotaci�n angular 3.805(Lobitz),2.114*2.14(Lobitz2005) 
TSR=7 #Tip speed ratio
Uinf=wy*35/TSR #velocidad libre
alfa=90/180*pi #Angulo de ataque libre
%Ttotal=0.5*(2*pi/wy) #Tiempo total de simulacion en segundos
Ttotal=50;


qf=-tan(0/180*pi); %pos flap inicial
qfdot=0; %vel flap inicial
amplitudflap=0/180*pi;
kflap=0.5; %quasi-steady 0.02, unsteady 0.1 and highly unsteady 0.5
frecuenciaflap=(2*Uinf*kflap/1)*(2*pi);

################################################################################
#Modelo aerodinamico############################################################
pitchangle=5*pi/180; #[en radianes] definen positivo antihorario. (2.6 optimo WindPACT)
pala.aero=csvread('input/palaWINDPACT_aero_rev1.csv'); #Sta number, 
raeroimport=pala.aero(:,1);
cuerdaimport=pala.aero(:,2);
twistimport=pala.aero(:,3)*pi/180+pitchangle; #CUIDADO Los datos definen positivo antihorario.
zpitchimport=pala.aero(:,4);



#Discretizado sabana adherida##############################################
m=10; #Divisiones en la cuerda, considero 6 para geomrotorupdate_r2!!!
n=25; #Divisiones en la envergadura, considero 20 para geomrotorupdate_r2!!!
desprendimiento=1; #1: borde de fuga, 2: borde de fuga y puntera derecha, 3: borde de fuga y ambas punteras

espaciado=2; #Metodo de posicionado de estaciones en envergadura:0 uniforme, 1 refinamiento senoidal hacia la punta, 2 refinamiento senoidal, 4 pala WindPACT.
[raero,n]=distr_nodos(espaciado,raeroimport,n,m);
numelm=m*n #Cantidad total de elementos

%pitchaxisdata=csvread('input/palaWINDPACT_str_revisado.csv'); %[sta,xestr,EA,EJzz,EJyy,EJyz,GJp,BMass,Izz,Iyy,Iyz,Sy,Sz,...,pitch];
%zpitch=interp1(pitchaxisdata(:,2)+h,pitchaxisdata(:,15),raero); #columna 2 es radio y columna 15 es fraccion de cuerda desde borde ataque
%clear pitchaxisdata;

cuerda=interp1(raeroimport,cuerdaimport,raero);
twist=interp1(raeroimport,twistimport,raero);
zpitch=interp1(raeroimport,zpitchimport,raero);

airfoilimport=csvread('input/airfoilnumber_WindPACT.csv');
airfoilnumber=interp1(airfoilimport(:,1),airfoilimport(:,2),raero,'nearest');

af1=csvread('input/S818_media.csv');
af2=csvread('input/S825_media.csv');
af3=csvread('input/S826_media.csv');

airfoildata=[af1,af2(:,2),af3(:,2)];#z vs y1 vs y2 vs y3

#Adimensionalizacion
Supaero=trapz(raeroimport,cuerdaimport)
mac=1/Supaero*trapz(raeroimport,cuerdaimport.^2) #cuerda media aerodinamica http://en.wikipedia.org/wiki/Chord_(aircraft)
LC=mac/10;	#Longitud caracteristica
xmac = interp1(cuerdaimport,raeroimport,mac)
xc=xmac;
%xc=0.75*(L+h); #maxima presicion en seccion tipica
VC=sqrt( (wy*(1.75+xc))^2 + (2/3*Uinf)^2 ) #velocidad en caracteristica
TC=LC/VC #Tiempo caracteristico
wyadim=wy*TC; #wy adim = wy dim  [rad/s]* TC[s]    wy dim = wy adim[rad]/TC[seg]
VlibreN=(Uinf/VC)*[0,sin(alfa),cos(alfa)]; #Adimensionalizado por VC
tsteps=floor(Ttotal/TC) #Pasos de tiempo de simulacion

#Creo malla y asigno posicion en el sistema B
[R0,elem]=geomrotorcreate_r4(m,n,raero,cuerda,twist,zpitch,airfoilnumber,airfoildata);

################################################################################
#Modelo estructural
h=1.75;#radio hub
L=35-h;#L pala elastica
#Inicializo modulo estructural
[M,Kc,Fc,Ipolar,Ke,Ktg,KT,nax,nflap,ntor,qe,autoval,autovec]=struct_init_r6(wy,L,h);

%KT=Ke; %PRUEBA!!!!!!!

Minv=inv(M);
LAMBDA=TC^2*Minv*KT; #Frecuencias reducidas por la velocidad caracteristica (TC=LC/VC)
A=[zeros(2*nflap+ntor),eye(2*nflap+ntor);-LAMBDA,zeros(2*nflap+ntor)];

#Condiciones iniciales para perturbaciones de coord generalizadas adimensionales en el tiempo
%Z0=[qe;zeros(2*nflap+nax,1)*TC];
Z0=[zeros(2*nflap+nax,1);zeros(2*nflap+nax,1)*TC];
Z0dot=[zeros(2*nflap+nax,1)*TC;zeros(2*nflap+nax,1)*(TC^2)];

%#Energia
%q=Z0(1:2*nflap+nax,1);
%qdot=Z0(2*nflap+nax+1:end,1)/TC;
%T=1/2*qdot'*M*qdot;
%U=1/2*q'*KT*q;
%Wnc(1)=0;
%ET(1)=T+U;
 
########################################################################################
#Inicializo variables
azimuth=0;
nodoswakeN=zeros(n,3);
nodoswakeB=zeros(n,3);
conwake=[ones(n,4),zeros(n,1)];
G=zeros(numelm,1);	
convect=0;
arranque=1;
initestela=0;
thetaxform=[0;azimuth;0]; #�ngulo para transformar entre sistemas coordenados

########################################################################################
#Posiciono malla en estado inicial
%[nodosB,elem]=geomrotorupdate_r4(n,m,R0,elem,[Z0(1:2);Z0(3:4)/TC],nax,nflap,ntor,wy,LC,TC,L,qf,qfdot,h);#Z debe ser dimensional
%[nodosB,elem]=geomrotorupdate_r4(n,m,R0,elem,[Z0(1:2*nflap+nax,1);Z0(2*nflap+nax+1:end,1)/TC],nax,nflap,ntor,wy,LC,TC,L,qf,qfdot,h);#Z debe ser dimensional
[nodosB,elem]=geomrotorupdate_r4(n,m,R0,elem,[qe+Z0(1:2*nflap+nax,1);Z0(2*nflap+nax+1:end,1)/TC],nax,nflap,ntor,wy,LC,TC,L,qf,qfdot,h);#Z debe ser dimensional
#Guardo malla inicial
nombre=strcat('output/blademesh.msh');
gmsh_export(nodosB*LC,elem,nodoswakeB*LC,conwake,G,0,m,n,nombre,0);
#############  
#Determino G inicial en tiempo 0 (sin estela) #No hay G para calcular la derivada y no hay estela
VlibreB=xformNB(-thetaxform,VlibreN);
G_old=G; %guardo valores viejos para calcular la variacion temporal de G
Ri=ones(1,3);
Rj=ones(1,3);
Rk=ones(1,3);
Rl=ones(1,3);
Gwake=0;
[G,Vindvecfull]=SOL_PHI_vec4(nodosB,elem,VlibreB,Ri,Rj,Rk,Rl,Gwake,numelm,arranque); %calculo circulaciones
DGDt=zeros(size(G,1),1);
#Evaluación de coeficientes y cargas
[dataout]=pressure_sp11(nodosB,elem,numelm,VlibreB,G,DGDt,nodoswakeB,conwake,m,n,Vindvecfull,desprendimiento); 
[Faero(:,1),FY3D(1),FY2D(1,:),FZ3D(1),FZ2D(1,:),MX3D(1),MX2D(1,:)]=cargasaero9(VC,LC,dataout,elem,nodosB,mac,Supaero,Uinf,L,rho,wy,nax,nflap,n,m);

arranque=0; #Ya arranco la estela, puedo calcular dG/dt
initestela=1; #A�n no convecte la estela
convect=1; #Convectar en el proximo paso
intervalo=1:tsteps;
elseif restart==1
  load aeroinit.mat;
  tsteps=length(FY3D)-1
  intervalo=tsteps:tsteps+500;
endif

% SOLO PARA PRUEBA ELASTICA
Faerop=zeros(2*nflap+nax,1);%prueba modelo elastico
Faero=zeros(2*nflap+nax,intervalo(end));
%Z0(1:2*nflap+nax,1)=0.1*autovec(:,14); #si nax1:modo10, nax2:modo13, nax3:modo14
Z0(1,1)=1;
% SOLO PARA PRUEBA ELASTICA

##########################################################################################
##Inicio simulacion
deltaTold=deltaT;
tiempo(1,1)=0;

for t=intervalo; #t es el tiempo siguiente.
fprintf("PASO DE TIEMPO %i\n",t)

%tarranque=-40; %10 cuerdas
%if t<tarranque
%    Faero=Faero*0;
%    modcarga=0;
%  elseif t<tarranque+40
%    modcarga=(t-tarranque)/40;
%  else
%    modcarga=1;
%endif
modcarga=1;


 ##################Fase Prediccion de nuevo estado###########################################################  
 #Avanzo estela a t+dt
%	azimuth=azimuth+wyadim*deltaT;
%	thetaxform=[0;azimuth;0]; #�ngulo para transformar entre marcos de referencia B y N
%	VlibreB=xformNB(-thetaxform,VlibreN);
%  [nodoswakeB,conwake]=convectar_rev14(nodosB,elem,G,nodoswakeB,conwake,VlibreB,deltaT,m,n,initestela,wyadim,desprendimiento);
%  initestela=0; #Ya convect� la estela por primera vez
%  G_old=G; %guardo circulacion anterior para calcular DGDt

%  #Avanzo consigna control
%  [qfnueva,qfdotnueva]=controlflap(t*TC,TC,Uinf);
%  qfplot(t+1)=-qfnueva*180/pi; #Guardo para graficar

  ##############################################################################
  ###############Interaccion fluido-estructura#####################################################

  #Evaluo aerodinamica en t+dt
      #Actualización de geometría sábana adherida 
      #Evaluación de circulaciones 
      #Evaluación de cargas
%[Faerop,FY3Dp,FY2Dp,FZ3Dp,FZ2Dp,MX3Dp,MX2Dp,nodosB,elem,nodoswakeB,G,dataout]=evalfun(n,m,R0,elem,[qe;zeros(2*nflap+nax,1)]+Z0(:,t),nax,nflap,ntor,wy,LC,TC,L,qfnueva,qfdotnueva,h,desprendimiento,nodoswakeB,conwake,VlibreB,numelm,arranque,G_old,deltaT,deltaTold,nodosB,VC,mac,Supaero,Uinf,rho);
%  Faerop=Faerop*modcarga;
  
  %Avanzo estructura a t+dt
%  [Zp]=struct_newmark(Z0(1:2*nflap+nax,t),Z0(2*nflap+nax+1:end,t),M,KT,Faerop,deltaT,TC);
  [Zp]=struct_HHT_Rev1(Z0(1:2*nflap+nax,t),Z0(2*nflap+nax+1:end,t),M,KT,Faerop,Faero(:,t),deltaT,TC);
 
  #Evaluo aerodinamica en t+dt
      #Actualización de geometría sábana adherida 
      #Evaluación de circulaciones 
      #Evaluación de cargas
%  [Faerop,FY3Dp,FY2Dp,FZ3Dp,FZ2Dp,MX3Dp,MX2Dp,nodosB,elem,nodoswakeB,G,dataout]=evalfun(n,m,R0,elem,[qe;zeros(2*nflap+nax,1)]+Zp,nax,nflap,ntor,wy,LC,TC,L,qfnueva,qfdotnueva,h,desprendimiento,nodoswakeB,conwake,VlibreB,numelm,arranque,G_old,deltaT,deltaTold,nodosB,VC,mac,Supaero,Uinf,rho);
%  Faerop=Faerop*modcarga;
  
 
  ##################Fase Correccion###########
  errZ=1;
  erraero=1;
  iter=0;
  while (errZ>1E-12)

    #Corrijo estructura con Newmark
%    [Zc]=struct_newmark(Z0(1:2*nflap+nax,t),Z0(2*nflap+nax+1:end,t),M,KT,Faerop,deltaT,TC);
    [Zc]=struct_HHT_Rev1(Z0(1:2*nflap+nax,t),Z0(2*nflap+nax+1:end,t),M,KT,Faerop,Faero(:,t),deltaT,TC);
    
    #Evaluo aerodinamica en t+dt
      #Actualización de geometría sábana adherida 
      #Evaluación de circulaciones 
      #Evaluación de cargas
%    [Faeroc,FY3Dp,FY2Dp,FZ3Dp,FZ2Dp,MX3Dp,MX2Dp,nodosB,elem,nodoswakeB,G,dataout]=evalfun(n,m,R0,elem,[qe;zeros(2*nflap+nax,1)]+Zc,nax,nflap,ntor,wy,LC,TC,L,qfnueva,qfdotnueva,h,desprendimiento,nodoswakeB,conwake,VlibreB,numelm,arranque,G_old,deltaT,deltaTold,nodosB,VC,mac,Supaero,Uinf,rho); %Evaluacion aero en t+dt
%    Faeroc=Faeroc*modcarga;
  
    errZ=norm(Zc-Zp,inf); %error estructural
%    erraero=norm((Faeroc-Faerop),inf)
    
    #Actualizacion para próxima iteracion 
%    Faerop=Faeroc;    
    Zp=Zc;
    iter++;
    
  endwhile
  iteraciones(t+1,1)=iter;
%  errorcarga(t+1,1)=erraero;
  errorconfiguracion(t+1,1)=errZ;
  ##################Fin avance DT###########################################################  

% Ajuste automatico del paso de tiempo 
%  deltaTold=deltaT;
%  if iter>20
%    deltaT=deltaT/2;
%    TC=TC/2;
%  elseif iter<10
%    deltaT=2*deltaT;
%    TC=TC*2;
%  endif
%  deltaT
%  TC
  tiempo(t+1,1)=tiempo(t,1)+deltaT*TC;
  
  #Actualizacion
  Z0dot(:,t+1)=[zeros(2*nflap+ntor,1);(TC^2)*Minv*Faerop]+A*Zp;
  Z0(:,t+1)=Zp;
%  Faero(:,t+1)=Faerop;
%  FY3D(t+1)=FY3Dp;
%  FZ3D(t+1)=FZ3Dp;
%  MX3D(t+1)=MX3Dp;
%  FY2D(t+1,:)=FY2Dp;
%  FZ2D(t+1,:)=FZ2Dp;

  
%  Pdin=0.5*rho*VC^2;
%  Aref=Supaero;
%  CN(t+1)=FY3D(t+1)/(Pdin*Aref); #La fuerza es la total sobre toda la superficie y la presion dinamica es [N]/[m2]
%  CDi(t+1)=FZ3D(t+1)/(Pdin*Aref); #La fuerza es la total sobre toda la superficie y la presion dinamica es [N]/[m2]
%  CM(t+1)=MX3D(t+1)/(Pdin*Aref*mac);
%  Lift(t+1)=FY3D(t+1)*cos(alfa)-FZ3D(t+1)*sin(alfa);
%  Drag(t+1)=FY3D(t+1)*sin(alfa)+FZ3D(t+1)*cos(alfa);
%  CL(t+1)=Lift(t+1)/(Pdin*Aref);
%  CDi(t+1)=Drag(t+1)/(Pdin*Aref);
  
   
%  #Salida de resultados###########################################################  
%	#Exporto coordenadas nodos en sistema B
%	nombre=strcat('output/step',int2str(t),'.msh');
%	gmsh_export(nodosB*LC,elem,nodoswakeB*LC,conwake,G,0,m,n,nombre,t)
%	nombre=strcat('output/bladecp',int2str(t),'.msh');
%	gmsh_export_bladecp(nodosB*LC,elem,G,dataout,m,n,nombre,t)
%%	[CNfinal,CMc4final,Xmomfinal,Ymomfinal,xsta,FYstafinal,FZstafinal,Empujefinal,Torquefinal,Potenciafinal]=postproceso_aero (VC,LC,dataout,elem,nodosB,mac,Supaero,Uinf,L,rho,wy,n,m,nodoswakeB,conwake,VlibreB,G,alfa);
%  nombre=strcat('output/ffsta_m',int2str(m),'_n',int2str(n),'_t',int2str(t),'.csv');
%  fid = fopen (nombre,'w');
%  xsta=[elem(1:n,5)*LC]'; #xcp para primera hilera de anillos
%  dlmwrite(fid,[xsta;FY2D(end,:);FZ2D(end,:)],' '); #Guardo distribucion de cargas en t final
%  fclose (fid);
  ################################################################################
  endfor



#FIN SIMULACION########################################################
#######################################################################



