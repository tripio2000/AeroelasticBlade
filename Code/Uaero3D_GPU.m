## Modelo aeroelastico
## Author: NGT

##Rev0: 9/2/19 Unifico módulos. Uso de plantilla aero3D_ala_rev4
##Rev1: 9/2/19 Agrego subiteraciones en estructura para trabajar con funciones de forma de alta frecuencia
##Rev2: ./2/19 
##Rev3: 20/2/19 Uso integrador Newmark-beta para estructura en vez de subiteraciones (m�s r�pido y estable)
##Rev4: 6/4/19 Agrego modulo aleron
##Rev5: 14/6/19 Uso estructura para datos de pala, aleron, etc
##Rev6: 18/3/22 Utilizo GPU para calcular velocidades inducidas

warning ("off", "Octave:broadcast");
warning ("off", "Octave:possible-matlab-short-circuit-operator");

clear all
close all
clc

pkg load ocl #paquete que habilita operaciones en GPU


restart=1;

if restart==0
###########################################################
#Parametros para integradores temporales
deltaT=1; #Siempre 1
#Parametros simulacion 
rho=1.225; #densidad del fluido, aire SI
wy=1*2.14 #velocidad rotaci�n angular 3.805(Lobitz),2.114*2.14(Lobitz2005) 
TSR=7 #Tip speed ratio
Uinf=wy*35/TSR #velocidad libre
alfa=90/180*pi #Angulo de ataque libre
%Ttotal=0.5*(2*pi/wy) #Tiempo total de simulacion en segundos
Ttotal=1.2;

tic

################################################################################
#Modelo aerodinamico############################################################
pala.pitchangle=5*pi/180; #[en radianes] definen positivo antihorario. (2.6 optimo WindPACT)
pala.aero=csvread('input/palaWINDPACT_aero_rev1.csv'); #Sta number, 

pala.Rtip=pala.aero(end,1); %radio rotor
pala.h=1.75;
pala.L=pala.Rtip-pala.h;

#Discretizado sabana adherida##############################################
pala.m=20; #Divisiones en la cuerda, considero 6 para geomrotorupdate_r2!!!
pala.n=30; #Divisiones en la envergadura, considero 20 para geomrotorupdate_r2!!!
desprendimiento=1; #1: borde de fuga, 2: borde de fuga y puntera derecha, 3: borde de fuga y ambas punteras

espaciado=2; #Metodo de posicionado de estaciones en envergadura:0 uniforme, 1 refinamiento senoidal hacia la punta, 2 refinamiento senoidal, 4 pala WindPACT.
[pala.raero,pala.n]=AeroDistrNodos(espaciado,pala.aero(:,1),pala.n,pala.m);
pala.numelm=pala.m*pala.n; #Cantidad total de elementos

pala.cuerda=interp1(pala.aero(:,1),pala.aero(:,2),pala.raero);
pala.twist=interp1(pala.aero(:,1),pala.aero(:,3)*pi/180+pala.pitchangle,pala.raero);
pala.zpitch=interp1(pala.aero(:,1),pala.aero(:,4),pala.raero);

airfoilimport=csvread('input/airfoilnumber_WindPACT.csv');
pala.airfoilnumber=interp1(airfoilimport(:,1),airfoilimport(:,2),pala.raero,'nearest');

af1=csvread('input/S818_media.csv');
af2=csvread('input/S825_media.csv');
af3=csvread('input/S826_media.csv');

pala.airfoildata=[af1,af2(:,2),af3(:,2)];#z vs y1 vs y2 vs y3

#Adimensionalizacion
pala.Supaero=trapz(pala.aero(:,1),pala.aero(:,2));
pala.mac=1/pala.Supaero*trapz(pala.aero(:,1),pala.aero(:,2).^2); #cuerda media aerodinamica http://en.wikipedia.org/wiki/Chord_(aircraft)
LC=pala.mac/5;	#Longitud caracteristica
xmac = interp1(pala.aero(:,2),pala.aero(:,1),pala.mac);
xc=xmac;
%xc=0.75*(L+h); #maxima presicion en seccion tipica
VC=sqrt( (wy*(1.75+xc))^2 + (2/3*Uinf)^2 ); #velocidad en caracteristica
TC=LC/VC #Tiempo caracteristico
wyadim=wy*TC; #wy adim = wy dim  [rad/s]* TC[s]    wy dim = wy adim[rad]/TC[seg]
VlibreN=(Uinf/VC)*[0,sin(alfa),cos(alfa)]; #Adimensionalizado por VC
tsteps=floor(Ttotal/TC) #Pasos de tiempo de simulacion

#Creo malla y asigno posicion en el sistema B
[R0,elem]=AeroGeomRotorCreate_r6(pala);

#################Modulo estructural#######################################################
##h=1.75;#radio hub
##L=35-h;#L pala elastica
#Inicializo modulo estructural
[M,Kc,Fc,Ipolar,Ke,Ktg,KT,nax,nflap,ntor,qe,autoval,autovec]=struct_init_r6(wy,pala.L,pala.h);

%KT=Ke; %PRUEBA!!!!!!!

Minv=inv(M);
LAMBDA=TC^2*Minv*KT; #Frecuencias reducidas por la velocidad caracteristica (TC=LC/VC)
A=[zeros(2*nflap+ntor),eye(2*nflap+ntor);-LAMBDA,zeros(2*nflap+ntor)];

#Condiciones iniciales para perturbaciones de coord generalizadas adimensionales en el tiempo
%Z0=[qe;zeros(2*nflap+nax,1)*TC];
Z0=[zeros(2*nflap+nax,1);zeros(2*nflap+nax,1)*TC];

%Z0(1,1)=0.1;

Z0dot=[zeros(2*nflap+nax,1)*TC;zeros(2*nflap+nax,1)*(TC^2)];

########################################################################################
#Inicializo variables
azimuth=0;
nodoswakeN=zeros(pala.n,3);
nodoswakeB=zeros(pala.n,3);
conwake=[ones(pala.n,4),zeros(pala.n,1)];
G=zeros(pala.numelm,1);	
convect=0;
arranque=1;
initestela=0;
thetaxform=[0;azimuth;0]; #�ngulo para transformar entre sistemas coordenados

########################################################################################
#Posiciono malla en estado inicial
[nodosB,elem]=AeroGeomRotorUpdate_r6(pala.n,pala.m,R0,elem,[qe+Z0(1:2*nflap+nax,1);Z0(2*nflap+nax+1:end,1)/TC],nax,nflap,ntor,wy,LC,TC,pala.L);#Z debe ser dimensional

#Guardo malla inicial
nombre=strcat('output/blademesh.msh');
gmsh_export(nodosB*LC,elem,nodoswakeB*LC,conwake,G,0,pala.m,pala.n,nombre,0);
#############  
#Determino G inicial en tiempo 0 (sin estela) #No hay G para calcular la derivada y no hay estela
VlibreB=xformNB(-thetaxform,VlibreN);
G_old=G; %guardo valores viejos para calcular la variacion temporal de G
Ri=ones(1,3);
Rj=ones(1,3);
Rk=ones(1,3);
Rl=ones(1,3);
Gwake=0;
[G,Vindvecfull]=AeroSolPHI_rev4(nodosB,elem,VlibreB,Ri,Rj,Rk,Rl,Gwake,pala.numelm,arranque); %calculo circulaciones
DGDt=zeros(size(G,1),1);
#Evaluación de coeficientes y cargas
[dataout]=AeroPressureSP_rev11(nodosB,elem,pala.numelm,VlibreB,G,DGDt,nodoswakeB,conwake,pala.m,pala.n,Vindvecfull,desprendimiento); 
[Faero(:,1),FY3D(1),FY2D(1,:),FZ3D(1),FZ2D(1,:),MX3D(1),MX2D(1,:)]=AeroLoads_rev10(VC,LC,dataout,elem,nodosB,pala.mac,pala.Supaero,Uinf,pala.L,rho,wy,nax,nflap,pala.n,pala.m);
 
arranque=0; #Ya arranco la estela, puedo calcular dG/dt
initestela=1; #A�n no convecte la estela
convect=1; #Convectar en el proximo paso
intervalo=1:tsteps;


elseif restart==1
  load aeroinit.mat;
  tsteps=length(FY3D)-1
  intervalo=tsteps:tsteps+1000;
endif

##########################################################################################
##Inicio simulacion
deltaTold=deltaT;
tiempo(1,1)=0;

for t=intervalo; #t es el tiempo siguiente.
  fprintf("PASO DE TIEMPO %i\n",t)

 ##################Fase Prediccion de nuevo estado###########################################################  
 ############################################################################################################  
 #Avanzo estela a t+dt#############################################################################
	azimuth=azimuth+wyadim*deltaT;
	thetaxform=[0;azimuth;0]; #�ngulo para transformar entre marcos de referencia B y N
	VlibreB=xformNB(-thetaxform,VlibreN);

  [nodoswakeB,conwake]=AeroConvect_rev14(nodosB,elem,G,nodoswakeB,conwake,VlibreB,deltaT,pala.m,pala.n,initestela,wyadim,desprendimiento);
  
  initestela=0; #Ya convect� la estela por primera vez
  G_old=G; %guardo circulacion anterior para calcular DGDt

  q=Z0(1:2*nflap+nax,t);
  dqdt=Z0(2*nflap+nax+1:end,t)/TC;
  d2qdt2=Z0dot(2*nflap+nax+1:end,t)/(TC^2);
  rhoaire=rho;
  ###############Interaccion fluido-estructura#####################################################

  #Evaluo aerodinamica en t+dt: Actualización de geometría sábana adherida , Evaluación de circulaciones , Evaluación de cargas
 [Faerop,FY3Dp,FY2Dp,FZ3Dp,FZ2Dp,MX3Dp,MX2Dp,nodosB,elem,nodoswakeB,G,dataout]=AeroEvalFun_r1(dataout,pala,R0,elem,[qe;zeros(2*nflap+nax,1)]+Z0(:,t),Z0dot(2*nflap+nax+1:end,t),nax,nflap,ntor,wy,LC,TC,desprendimiento,nodoswakeB,conwake,VlibreB,arranque,G_old,deltaT,deltaTold,nodosB,VC,Uinf,rho);

  
  %Avanzo estructura a t+dt
  [Zp,d2qdt2]=struct_HHT_Rev1(Z0(1:2*nflap+nax,t),Z0(2*nflap+nax+1:end,t),M,KT,Faerop,Faero(:,t),deltaT,TC);

 
  #Evaluo aerodinamica en t+dt: Actualización de geometría sábana adherida , Evaluación de circulaciones , Evaluación de cargas
 [Faerop,FY3Dp,FY2Dp,FZ3Dp,FZ2Dp,MX3Dp,MX2Dp,nodosB,elem,nodoswakeB,G,dataout]=AeroEvalFun_r1(dataout,pala,R0,elem,[qe;zeros(2*nflap+nax,1)]+Zp,d2qdt2,nax,nflap,ntor,wy,LC,TC,desprendimiento,nodoswakeB,conwake,VlibreB,arranque,G_old,deltaT,deltaTold,nodosB,VC,Uinf,rho);

  ##################Fase Correccion##############################
  errZ=1;
  erraero=1;
  iter=0;
  
  while (errZ>1E-12)&&iter<10

    #Corrijo estructura con Newmark
    [Zc,d2qdt2]=struct_HHT_Rev1(Z0(1:2*nflap+nax,t),Z0(2*nflap+nax+1:end,t),M,KT,Faerop,Faero(:,t),deltaT,TC);
    
    #Evaluo aerodinamica en t+dt: Actualización de geometría sábana adherida , Evaluación de circulaciones , Evaluación de cargas
   [Faeroc,FY3Dp,FY2Dp,FZ3Dp,FZ2Dp,MX3Dp,MX2Dp,nodosB,elem,nodoswakeB,G,dataout]=AeroEvalFun_r1(dataout,pala,R0,elem,[qe;zeros(2*nflap+nax,1)]+Zc,d2qdt2,nax,nflap,ntor,wy,LC,TC,desprendimiento,nodoswakeB,conwake,VlibreB,arranque,G_old,deltaT,deltaTold,nodosB,VC,Uinf,rho); %Evaluacion aero en t+dt

  
    errZ=norm(Zc-Zp,inf); %error estructural
    erraero=norm((Faeroc-Faerop),inf)
    
    #Actualizacion para próxima iteracion 
    Faerop=Faeroc;    
    Zp=Zc;
    iter++;
  endwhile
  iteraciones(t+1,1)=iter;
  errorcarga(t+1,1)=erraero;
  errorconfiguracion(t+1,1)=errZ;
  #Actualizacion  
  tiempo(t+1,1)=tiempo(t,1)+deltaT*TC;
%  Z0dot(:,t+1)=[zeros(2*nflap+ntor,1);(TC^2)*Minv*Faerop]+A*Zp;
  Z0dot(:,t+1)=[Zp(2*nflap+ntor+1:end);d2qdt2];
  Z0(:,t+1)=Zp;
  Faero(:,t+1)=Faerop;
  ##################Fin avance DT#############################################################################   
  ############################################################################################################     

  #Salida de resultados###########################################################  
  FY3D(t+1)=FY3Dp;
  FZ3D(t+1)=FZ3Dp;
  MX3D(t+1)=MX3Dp;
  FY2D(t+1,:)=FY2Dp;
  FZ2D(t+1,:)=FZ2Dp;
  Pdin=0.5*rho*VC^2;
  Aref=pala.Supaero;
  CN(t+1)=FY3D(t+1)/(Pdin*Aref); #La fuerza es la total sobre toda la superficie y la presion dinamica es [N]/[m2]
  CDi(t+1)=FZ3D(t+1)/(Pdin*Aref); #La fuerza es la total sobre toda la superficie y la presion dinamica es [N]/[m2]
  CM(t+1)=MX3D(t+1)/(Pdin*Aref*pala.mac);
  Lift(t+1)=FY3D(t+1)*cos(alfa)-FZ3D(t+1)*sin(alfa);
  Drag(t+1)=FY3D(t+1)*sin(alfa)+FZ3D(t+1)*cos(alfa);
  CL(t+1)=Lift(t+1)/(Pdin*Aref);
  CDi(t+1)=Drag(t+1)/(Pdin*Aref);
	#Exporto coordenadas nodos en sistema B
	nombre=strcat('output/step',int2str(t),'.msh');
	gmsh_export(nodosB*LC,elem,nodoswakeB*LC,conwake,G,0,pala.m,pala.n,nombre,t)
	nombre=strcat('output/bladecp',int2str(t),'.msh');
	gmsh_export_bladecp(nodosB*LC,elem,G,dataout,pala.m,pala.n,nombre,t)
%	[CNfinal,CMc4final,Xmomfinal,Ymomfinal,xsta,FYstafinal,FZstafinal,Empujefinal,Torquefinal,Potenciafinal]=postproceso_aero (VC,LC,dataout,elem,nodosB,mac,Supaero,Uinf,L,rho,wy,n,m,nodoswakeB,conwake,VlibreB,G,alfa);
  nombre=strcat('output/ffsta_m',int2str(pala.m),'_n',int2str(pala.n),'_t',int2str(t),'.csv');
  fid = fopen (nombre,'w');
  xsta=[elem(1:pala.n,5)*LC]'; #xcp para primera hilera de anillos
  dlmwrite(fid,[xsta;FY2D(end,:);FZ2D(end,:)],' '); #Guardo distribucion de cargas en t final
  fclose (fid);

  endfor

#FIN SIMULACION########################################################
#######################################################################



