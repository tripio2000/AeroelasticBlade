%postproceso FSI WindPACT
clc
close all


#puntos a graficar
xb=linspace(h,L+h,20);
#evaluo funciones de forma
[ui,uip,uipp,vi,vip,vipp,wi,wip,wipp,ti,tip,tipp,nax,nflap,nedge,ntor]=modoseval_r1(xb,L,nax,nflap);
%las funciones de forma se almacenan como: [modos x posicion axial]

#coord generalizadas adimensionales en el tiempo
qv=Z0(1:nflap,:);
qw=Z0(nflap+1:2*nflap,:);
qtx=Z0(2*nflap+1:2*nflap+nax,:);

#desplazamientos del centro elástico [posicion x tiempo]
v0=vi'*qv; #vi * qv
w0=wi'*qw;#wi * qw
thetax=ti'*qtx;#ti * qt
thetay=vip'*qv;#ti * qt
thetaz=wip'*qw;#ti * qt

#Velocidades
qend=2*nflap+nax;
dotqv=Z0(qend+1:qend+nflap,:);
dotqw=Z0(qend+nflap+1:qend+2*nflap,:);
dotqtx=Z0(qend+2*nflap+1:end,:);
	
dotv0=vi'*dotqv; #vi * dot qv
dotw0=wi'*dotqw; #wi * dot qw
dotthetax=ti'*dotqtx; #ti * dot qt
dotthetay=vip'*dotqv; #thetay=dv0/dx, vip dot qv
dotthetaz=wip'*dotqw; #thetaz=dw0/dx, #wip dot qw 
  
%figure(1)
%subplot(1,3,1)
%plot(v0(end,:),dotv0(end,:),'o-r','linewidth',2)
%xlabel('Desplazamiento de punta - v')
%ylabel('Velocidad')
%subplot(1,3,2)
%plot(w0(end,:),dotw0(end,:),'o-b','linewidth',2)
%xlabel('Desplazamiento de punta - w')
%ylabel('Velocidad')
%subplot(1,3,3)
%plot(thetax(end,:),dotthetax(end,:),'o-k','linewidth',2)
%xlabel('Desplazamiento de punta - tx')
%ylabel('Velocidad')
%
%figure(2)
%subplot(1,3,1)
%plot(tiempo,v0(end,:),'o-r','linewidth',2)
%xlabel('Paso de tiempo')
%ylabel('Desplazamiento de punta - v [m]')
%
%subplot(1,3,2)
%plot(tiempo,w0(end,:),'o-b','linewidth',2)
%xlabel('Paso de tiempo')
%ylabel('Desplazamiento de punta - w [m]')
%
%subplot(1,3,3)
%plot(tiempo,thetax(end,:)*180/pi,'o-k','linewidth',2)
%xlabel('Paso de tiempo')
%ylabel('Desplazamiento de punta - tx [grados]')
%
%puntos=length(CL);
%frecuencias = (0:puntos-1)/(puntos*TC); #en Hz
%amplitudesv=abs(fft(v0(1:puntos)));
%amplitudest=abs(fft(thetax(1:puntos)));
%amplitudesqv=abs(fft(qv')); %cada columna es una fft distinta
%amplitudesqtx=abs(fft(qtx')); %cada columna es una fft distinta


set (0, "defaultaxesfontname", "Times New Roman")
set (0, "defaulttextfontname", "Times New Roman") 
set (0, "defaultaxesfontsize", 11)
set (0, "defaulttextfontsize", 11)
graphics_toolkit('gnuplot') 

h=figure(1)
W = 5; H = 5;
ventanatiempo=1:columns(v0);
subplot(1,3,1)
plot(v0(end,ventanatiempo),dotv0(end,ventanatiempo),'-r','linewidth',1)
xlabel("Desplazamiento de punta\n v[m]")
ylabel("Velocidad\n dvdt[m]")
subplot(1,3,2)
plot(w0(end,ventanatiempo),dotw0(end,ventanatiempo),'-b','linewidth',1)
xlabel("Desplazamiento de punta\n w[m]")
ylabel('Velocidad - dwdt[m/s]')
subplot(1,3,3)
plot(thetax(end,ventanatiempo)*180/pi,dotthetax(end,ventanatiempo)*180/pi,'-k','linewidth',1)
xlabel('Desplazamiento de punta - tx[grados]')
ylabel('Velocidad - dtxdt [grados/s]')

set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print 'plano de fases.eps' -depsc -tight

h=figure(2)
W = 5; H = 7;
subplot(3,1,1)
##lambdav=0.006;%0.003(beta0.02) 0.012(beta0.1) 0.015(beta1/3)
##valormax=max(v0(end,ventanatiempo));
##tiempovalormax=find(v0(end,:)==valormax)
##plot(tiempo,v0(end,:),'-r','linewidth',2,(tiempovalormax*TC:tiempo(end)), valormax*e.^(-lambdav*((tiempovalormax*TC:tiempo(end))-tiempovalormax*TC)))
plot(tiempo,v0(end,:),'-r','linewidth',2)
##xlabel("Tiempo - t[s]")
ylabel("Desplazamiento de punta\n v[m]")
subplot(3,1,2)
##lambdaw=0.006;%0.003(beta0.02) 0.012(beta0.1) 0.015(beta1/3)
##valormax=max(w0(end,ventanatiempo));
##tiempovalormax=find(w0(end,:)==valormax);
##plot(tiempo,w0(end,:),'-b','linewidth',2,(tiempovalormax*TC:tiempo(end)), valormax*e.^(-lambdaw*((tiempovalormax*TC:tiempo(end))-tiempovalormax*TC)))
plot(tiempo,w0(end,:),'-b','linewidth',2)
##xlabel("Tiempo - t[s]")
ylabel("Desplazamiento de punta\n w[m]")
subplot(3,1,3)
##lambdat=0.014;%0.006(beta0.02) 0.014(beta0.1) 0.015(beta1/3)
##valormax=max(thetax(end,ventanatiempo));
##tiempovalormax=find(thetax(end,:)==valormax);
##plot(tiempo,thetax(end,:)*180/pi,'-k','linewidth',2,(tiempovalormax*TC:tiempo(end)),180/pi*valormax*e.^(-lambdat*((tiempovalormax*TC:tiempo(end))-tiempovalormax*TC)))
plot(tiempo,thetax(end,:)*180/pi,'-k','linewidth',2)
xlabel("Tiempo - t[s]")
ylabel("Desplazamiento de punta\n tx[grados]")

set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print 'evolucion.eps' -depsc -tight

h=figure(8)
W = 5; H = 7;
recorte=floor(1/TC) %muestro cambios en ultimo segundo

subplot(3,1,2)
plot(xb,w0(:,end-recorte:end),'linewidth',2)
##xlabel('posicion x_b[m]')
ylabel('desplazamiento - w[m]')
subplot(3,1,3)
plot(xb,thetax(:,end-recorte:end)*180/pi,'linewidth',2)
xlabel('posicion x_b[m]')
ylabel('Giro - tx[grados]')
subplot(3,1,1)
plot(xb,v0(:,end-recorte:end),'linewidth',2)
##xlabel('posicion x_b[m]')
ylabel('desplazamiento - v[m]')
title("Evolucion de la configuraci�n \n durante el ultimo segundo de simulacion")

set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print 'configuracion_ultimo_segundo.eps' -depsc2

%figure(9)
%subplot(1,3,1)
%plot(tiempo,Z0(1:nflap,:))
%title('Evolución de coordenadas generalizadas v')
%xlabel('Tiempo - t[s]')
%ylabel('Amplitud de la coordenada generalizada')
%subplot(1,3,2)
%plot(tiempo,Z0(nflap+1:2*nflap,:))
%title('Evolución de coordenadas generalizadas w')
%xlabel('Tiempo - t[s]')
%ylabel('Amplitud de la coordenada generalizada')
%subplot(1,3,3)
%plot(tiempo,Z0(2*nflap:end,:))
%title('Evolución de coordenadas generalizadas t')
%xlabel('Tiempo - t[s]')
%ylabel('Amplitud de la coordenada generalizada')




%figure(10)
%subplot(1,3,1)
%plot(tiempo,Faero(1:nflap,:))
%title('Evolución de cargas generalizadas v')
%xlabel('Paso de tiempo')
%ylabel('Amplitud de la carga generalizada')
%subplot(1,3,2)
%plot(tiempo,Faero(nflap+1:2*nflap,:))
%title('Evolución de cargas generalizadas w')
%xlabel('Paso de tiempo')
%ylabel('Amplitud de la carga generalizada')
%subplot(1,3,3)
%plot(tiempo,Faero(2*nflap:end,:))
%title('Evolución de cargas generalizadas t')
%xlabel('Paso de tiempo')
%ylabel('Amplitud de la carga generalizada')

puntos=columns(Z0);
tamventana=256;
i=0;
ventana=(puntos-(i+1)*(tamventana-1):puntos-i*(tamventana-1));
amplitudesv=abs(fft(v0(end,ventana)));
amplitudesw=abs(fft(w0(end,ventana)));
amplitudest=abs(fft(thetax(end,ventana)));
frecuencias=(0:(tamventana-1))/(tamventana*TC); #en Hz
fprintf('el df es %f[Hz]\n',1/(tamventana*TC))

h=figure(4)
W = 5; H = 5;
semilogy(frecuencias,amplitudesv,'linewidth',2,frecuencias,amplitudesw,'linewidth',2,frecuencias,amplitudest,'linewidth',2)
%%semilogy(frecuencias',amplitudesv,'x-','linewidth',2)
legend('v_o (tip)','w_o (tip)','theta_x (tip)')
xlabel("Frecuencia [Hz]")
ylabel("abs(fft(amplitud))")
axis([0,frecuencias(end)/2,min(amplitudest),max(amplitudesv)])

valorw=sort(amplitudesw(1:end/2));
frecfundw=frecuencias(amplitudesw==valorw(end-1))
valorv=max(amplitudesv(frecuencias>3&frecuencias<10));
frecfundv=frecuencias(amplitudesv==valorv)

text(frecfundw(1)+0.1,valorw(end-1),[num2str(frecfundw(1),3),'Hz'])
text(frecfundv(1)+0.1,valorv,[num2str(frecfundv(1),3),'Hz'])

grid on

set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print 'fft.eps' -depsc2

%figure(4)
%subplot(1,2,1)
%%semilogy(frecuencias,amplitudesv,'linewidth',2,frecuencias,amplitudest,'linewidth',2)
%semilogy(frecuencias',amplitudesv,'x-','linewidth',2)
%xlabel("Frecuencia [Hz]")
%ylabel("abs(fft(amplitudesv))")
%axis([0,30,0.01,1000])
%
%subplot(1,2,2)
%%semilogy(frecuencias,amplitudesv,'linewidth',2,frecuencias,amplitudest,'linewidth',2)
%semilogy(frecuencias',amplitudest,'x-','linewidth',2)
%xlabel("Frecuencia [Hz]")
%ylabel("abs(fft(amplitudest))")
%axis([0,30,0.01,1000])

xcp=elem(1:n,5)*LC;
Mz=dot(xcp,FY2D(end,:))
My=dot(xcp,FZ2D(end,:))
Torquelength=xcp.*FZ2D(end:end,:)';

#Datos WTPerf
%Relm  IncidAng Thr/Len	Trq/Len
dataWindPACT=csvread('WindPACT_WTPERF.csv');
%Elem	RElm	IncidAng	Azim	LocVel	Re	Loss	AxInd	TanInd	AFAngle	Alpha	Cl	Cd	Cm	ThrCo	TrqCo	PwrCo	Thr/Len	Trq/Len	Power


%figure(5)
%inflowangle=90-atan(FY2D(end,:)./abs(FZ2D(end,:)))*180/pi;
%
%
%plot(xb,thetax(:,end)*180/pi,'linewidth',2,xcp,inflowangle,'linewidth',2,raeroimport,twistimport*180/pi,'linewidth',2)
%xlabel('posicion x_b[m]')
%ylabel('ángulo [grados]')
%legend('giro elaśtico','inflowangle','alabeo aerodinamico')

%figure(6)
%plot(xcp,inflowangle,'linewidth',2,dataWindPACT(:,2),dataWindPACT(:,10),'linewidth',2)
%xlabel('posicion x_b[m]')
%ylabel('ángulo [grados]')
%legend('inflowangle UVLM','AFAngle WTPerf')
%
%figure(7)
%subplot(1,2,1)
%plot(xcp,FY2D(end:end,:),'linewidth',2,dataWindPACT(:,2),dataWindPACT(:,18),'linewidth',2)
%xlabel('posicion x_b[m]')
%ylabel('Thr/Len [N/m]')
%legend('UVLM','WTPerf','location','south')
%
%subplot(1,2,2)
%plot(xcp,-Torquelength,'linewidth',2,dataWindPACT(:,2),dataWindPACT(:,19),'linewidth',2)
%xlabel('posicion x_b[m]')
%ylabel('Torque/length [N]')
%legend('UVLM','WTPerf','location','south')

%
%figure(8)
%plot(xb,thetax(:,end:end)*180/pi,'linewidth',2)
%xlabel('posicion x_b[m]')
%ylabel('ángulo [grados]')


%	#########################################POSTPROCESO estructura######################################
%filename='palaWINDPACT_str_revisionfinal.csv'; %viga no uniforme con centros desplazados
%pala.estr=csvread(['input/',filename]); #Leo propiedades desde archivo externo
%
%#Organizo datos
%#Propiedades elasticas:
%xestr=pala.estr(:,2); 
%EA=pala.estr(:,3);
%EJzz=pala.estr(:,4);#flap
%EJyy=1E10*pala.estr(:,5);#edge multiplicar*1E10 para viga plana
%EJyz=pala.estr(:,6);
%GJp=pala.estr(:,7);% multiplicar*1E5 para viga plana
%#Propiedades inerciales:
%BMass=pala.estr(:,8);
%Izz=pala.estr(:,9);
%Iyy=pala.estr(:,10);
%Iyz=pala.estr(:,11);
%Sy=pala.estr(:,12); %si no considero acoplamiento inercial, multiplicar por 0
%Sz=0*pala.estr(:,13); %si no considero acoplamiento inercial, multiplicar por 0
%Itt=(Izz+ Iyy)';
%
%#Divido la viga en segmentos
%N=300; #Puntos de integracion del modelo estructural
%x=linspace(0,L,N);
%deltax=L/(N-1); #longitud de cada segmento a integrar
%xi=x(1,1:N-1); #punto inicial de cada segmento a integrar
%#Interpolo propiedades en los puntos de cada elemento
%xdata=xi;
%BMassinterp=interp1(xestr,BMass,xdata);
%EAinterp=interp1(xestr,EA,xdata);
%EJzzinterp=interp1(xestr,EJzz,xdata);
%EJyyinterp=interp1(xestr,EJyy,xdata);
%EJyzinterp=interp1(xestr,EJyz,xdata);
%GJpinterp=interp1(xestr,GJp,xdata);
%Izzinterp=interp1(xestr,Izz,xdata);
%Iyyinterp=interp1(xestr,Iyy,xdata);
%Iyzinterp=interp1(xestr,Iyz,xdata);
%Szinterp=interp1(xestr,Sz,xdata);
%Syinterp=interp1(xestr,Sy,xdata);
%#Evaluo funciones de forma en xdata
%puntosdata=size(xdata,2);
%[ui,uip,uipp,vi,vip,vipp,wi,wip,wipp,ti,tip,tipp,nax,nflap,nedge,ntor]=modoseval_r1(xdata,L,nax,nflap); %Evaluacion de funciones de forma fi
%  
%  %desplazamiento del centro elastico
%%	u0=qe(1:nax,1)'*ui; %q(t) es coord generalizada y ui es funcion de forma fi_u^i(x)
%%	v0=qe(1:nflap,1)'*vi;
%%	w0=qe(nflap+1:2*nflap,1)'*wi;
%%	theta0=qe(2*nflap+1:2end,1)'*ti;
%v0=vi'*qv; #vi * qv
%w0=wi'*qw;#wi * qw
%thetax=ti'*qtx;#ti * qt
%thetay=vip'*qv;#ti * qt
%thetaz=wip'*qw;#ti * qt
%  
%  %deformaciones
%%	varepsilon0=qe(1:nax,1)'*uip;
%%	thetay=-qe(nax+nflap+1:nax+2*nflap,1)'*wip;
%%	thetaz=qe(nax+1:nax+nflap,1)'*vip;
%%	psy=qe(nax+2*nflap+1:2*nax+2*nflap,1)'*tip;
%%	kappay=-qe(nax+nflap+1:nax+2*nflap,1)'*wipp;
%%	kappaz=qe(nax+1:nax+nflap,1)'*vipp;
%	psy=tip'*qtx;
%	kappay=-wipp'*qw;
%	kappaz=vipp'*qv;
%	
%  #Fuerzas seccionales @ x_estructural:
%%  axialstrain=varepsilon0+0.5*(thetay.^2+thetaz.^2);
%  axialstrain=0.5*(thetay(:,end).^2+thetaz(:,end).^2);
%	Nx=EAinterp.*axialstrain;
%	My=EJyyinterp.*kappay(:,end)-EJyzinterp.*kappaz(:,end); #Medge
%	Mz=EJyzinterp.*kappay(:,end)-EJzzinterp.*kappaz(:,end); #Mflap
%	Tx=GJpinterp.*psy(:,end);
 

