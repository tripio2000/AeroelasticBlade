function [M,Kc,Fc,Ipolar,Ke,Ktg,KT,nax,nflap,ntor,qe,autoval,autovec]=struct_init_r6(wy,L,h)
## Calculo de matrices del modelo estructural
## referencia: Modelo estructural de vigas 3D linealizado revisado3.odt y viga3D_rev18.m
## uso: [M,Kc,Fc,Ipolar,Ke,Ktg,nax,nflap,ntor,qe,xmedio,autoval,autovec]=struct_init_r4(filename,wy,L,h);
## filename='input/pala1.5MWpaperBir_rev1.csv'; %archivo con distribuciones de propiedades #Sta number,Zblade,EA,EJzz,EJyy,EJyz,GJp,Bmass,Izz,Iyy,Iyz,Sy,Sz
## wy=2.14 #velocidad rotaci�n angular [rad/s]
## L=33.25; #L pala [m]
## h=1.75; #radio hub [m]


## Author: NGT
## Rev1: Reemplazo funciones de forma por formas modales del libro de Hodges.
## Rev2: Depuracion. Elimino opcionales. Solo sirve para simular pala de tesis con funciones modales.
## Rev3: Leo parametros de aeroflex (reduzco errores de actualizacion)
## Rev4: Elimino dinamica de perturbaciones de grado axial
## Rev5: Reviso esquema de integracion
## Rev6: simplifico para usar en aeroflex, elimino postprocesos

%filename='pala1.5MWpaperBir_rev1.csv'; %viga no uniforme con centros desplazados
filename='palaWINDPACT_str_revisionfinal3.csv'; %viga no uniforme con centros desplazados

pala.estr=csvread(['input/',filename]); #Leo propiedades desde archivo externo

nax=1;
nflap=1;

nmodos=nax+2*nflap
pala.estr
#Organizo datos
#Propiedades elasticas:
xestr=pala.estr(:,2)-1.75; 
EA=pala.estr(:,3);
EJzz=pala.estr(:,4);#flap 0.5 para 1 modo
EJyy=pala.estr(:,5);#edge multiplicar*1E10 para viga plana
EJyz=pala.estr(:,6); #0.7 para 1 modo
GJp=pala.estr(:,7);% multiplicar*1E5 para viga plana, 0.25 para 1 modo
#Propiedades inerciales:
BMass=pala.estr(:,8);
Izz=pala.estr(:,9);
Iyy=pala.estr(:,10);
Iyz=pala.estr(:,11);
Sy=pala.estr(:,12); %si no considero acoplamiento inercial, multiplicar por 0
Sz=pala.estr(:,13); %si no considero acoplamiento inercial, multiplicar por 0
Itt=(Izz+Iyy)';

#Divido la viga en segmentos
N=300; #Puntos de integracion del modelo estructural
x=linspace(0,L,N);
deltax=L/(N-1); #longitud de cada segmento a integrar
xi=x(1,1:N-1); #punto inicial de cada segmento a integrar

#Interpolo propiedades en los puntos de Gauss de cada elemento
#primer punto de Gauss
x1=xi+0.211325*deltax*ones(1,N-1);
BMassinterp1=interp1(xestr,BMass,x1);
EAinterp1=interp1(xestr,EA,x1);
EJzzinterp1=interp1(xestr,EJzz,x1);
EJyyinterp1=interp1(xestr,EJyy,x1);
EJyzinterp1=interp1(xestr,EJyz,x1);
GJpinterp1=interp1(xestr,GJp,x1);
Izzinterp1=interp1(xestr,Izz,x1);
Iyyinterp1=interp1(xestr,Iyy,x1);
Iyzinterp1=interp1(xestr,Iyz,x1);
Szinterp1=interp1(xestr,Sz,x1);
Syinterp1=interp1(xestr,Sy,x1);
#Evaluo funciones de forma en punto de Gauss
[ui1,uip1,uipp1,vi1,vip1,vipp1,wi1,wip1,wipp1,ti1,tip1,tipp1,nax,nflap,nedge,ntor]=modoseval_r1(x1,L,nax,nflap);
#segundo punto de Gauss
x2=xi+0.788675*deltax*ones(1,N-1);
BMassinterp2=interp1(xestr,BMass,x2);
EAinterp2=interp1(xestr,EA,x2);
EJzzinterp2=interp1(xestr,EJzz,x2);
EJyyinterp2=interp1(xestr,EJyy,x2);
EJyzinterp2=interp1(xestr,EJyz,x2);
GJpinterp2=interp1(xestr,GJp,x2);
Izzinterp2=interp1(xestr,Izz,x2);
Iyyinterp2=interp1(xestr,Iyy,x2);
Iyzinterp2=interp1(xestr,Iyz,x2);
Szinterp2=interp1(xestr,Sz,x2);
Syinterp2=interp1(xestr,Sy,x2);
#Evaluo funciones de forma en punto de Gauss
[ui2,uip2,uipp2,vi2,vip2,vipp2,wi2,wip2,wipp2,ti2,tip2,tipp2,nax,nflap,nedge,ntor]=modoseval_r1(x2,L,nax,nflap);
%evaluo punto medio de cada segmento
xdata=0.5*x1+0.5*x2; 
[ui,uip,uipp,vi,vip,vipp,wi,wip,wipp,ti,tip,tipp,nax,nflap,nedge,ntor]=modoseval_r1(xdata,L,nax,nflap); %Evaluacion de funciones de forma fi


###########################################################################################################
#Evaluo las matrices lineales integrando con 2 puntos de gauss
###########################################################################################################
[M,Kc,Fc,Ipolar,Ke]=linmat2_r3(h,x1,x2,deltax,BMassinterp1,Izzinterp1,Iyyinterp1,Iyzinterp1,Syinterp1,Szinterp1,EAinterp1,EJzzinterp1,EJyzinterp1,EJyyinterp1,GJpinterp1,BMassinterp2,Izzinterp2,Iyyinterp2,Iyzinterp2,Syinterp2,Szinterp2,EAinterp2,EJzzinterp2,EJyzinterp2,EJyyinterp2,GJpinterp2,ui1,ui2,vi1,vi2,wi1,wi2,ti1,ti2,uip1,uip2,vip1,vip2,wip1,wip2,vipp1,vipp2,wipp1,wipp2,tip1,tip2,nax,nflap,nax);
#############################################################################################################
#Linealizacion de la matriz geométrica alrededor de estado de operacion
###########################################################################################################
#Determino punto de equilibrio con esquema incremental iterativo
%propongo la solucion lineal como punto inicial 
Klineal=(Ke - wy^2*Kc);
Fexterna=(wy^2*Fc);
qu=Klineal\Fexterna; 

#Busco punto de equilibrio por newton-raphson
for lambda=1:10 %incrementos de carga
	Fe=Fexterna*(lambda/10);
	errorDu=1;
	iter=1;
	while and(errorDu>1e-15,iter<50) %Solucion de sistema nolineal
		#Actualizo matriz tangente y esfuerzos internos
		[Ktg,Kg]=Kgeometrica3(x1,x2,qu,uip1,uip2,vip1,vip2,wip1,wip2,vipp1,vipp2,wipp1,wipp2,tip1,tip2,EAinterp1,EAinterp2,deltax,Ke,Kc,wy,nax,nflap,ntor);
		KT=(Ke - wy^2*Kc+Ktg); 
		Fi=(Ke - wy^2*Kc+Kg)*qu;
		#Encuentro incrementos
		Du = KT \ (Fe-Fi);
		errorDu=norm(Du,1);
		#Actualizo solucion
		qu=qu+Du;
	  iter=iter+1;	
	endwhile
	iter=iter-1;
endfor

qe=qu; %punto de equilibrio
#R4:Elimino grado de libertad axial (Desprecio dinamica de perturbacion axial frente a las transversales y torsional)
M(1:nax,:)=[];
M(:,1:nax)=[];
Kc(1:nax,:)=[];
Kc(:,1:nax)=[];
Ke(1:nax,:)=[];
Ke(:,1:nax)=[];
Ktg(1:nax,:)=[];
Ktg(:,1:nax)=[];
KT(1:nax,:)=[];
KT(:,1:nax)=[];
Fc=Fc(nax+1:end,1);
qe=qu(nax+1:end,1);

#Analisis modal alrededor del punto de equilibrio
[autovec,autoval]=eig(KT,M);
autovalvec=diag(autoval); %recupero vector de autovalores
autovalvecsort=sort(autovalvec); %ordeno autovalores de menor a mayor
frecgirorad=sqrt(diag(autoval));  #en [rad/s]
disp("Autovalores ordenados, expresados en rad/s")
frecordenadarad=sort(frecgirorad)
frecgiro=frecgirorad/(2*pi);  #en [Hz]
#Determino el paso de tiempo para simulaciones
frecordenada=sort(frecgiro);
tminimo=min(1./frecordenada); #Período más pequeño
tsim=tminimo/9;

for i=1:nmodos
  colautovector(i,:)=any(autoval==autovalvecsort(i)); %busco columna del modo i y guardo en colautovector
   #modos ordenados(de menor autovalor)
    autovectorordenado=autovec(:,colautovector(i,:));
%    modou0(i,:)=autovectorordenado(1:nax,1)'*ui; 
    modov0(i,:)=autovectorordenado(1:nflap,1)'*vi;
	  modow0(i,:)=autovectorordenado(nflap+1:2*nflap,1)'*wi;
	  modotheta0(i,:)=autovectorordenado(2*nflap+1:2*nflap+nax,1)'*ti;
endfor



endfunction

