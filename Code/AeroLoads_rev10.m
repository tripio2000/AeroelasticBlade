function [Faero,FY3D,FY2D,FZ3D,FZ2D,MX3D,MX2D]=AeroLoads_rev10(VC,LC,dataout,elem,nodosB,mac,Supaero,Uinf,L,rho,wy,nax,nflap,n,m)
## "ver Modelo estructural de vigas 3D linealizado revisado3.odt"
## Revisado calculo de Nv,Nw,Nt
## Vuelvo a 3D pala
## Elimino grado axial
## Rev6: reviso esquema de cuadratura usando aproximación polinomial para el campo de presión.
## Rev7: 10/1/19 reviso integrandos usando integrales de superficie
## Rev8: 9/2/19 cambio coord generalizadas para modelo WindPACT
## Rev9: 26/2/19 Integro cargas generalizadas con Gauss de 3 puntos.
## Rev10: 27/4/19 Mejoro calculo de fuerzas seccionales
#(VC,LC,dataout,elem,nodosB,poliu,poliv,polivprima,polit,Lflap,elemflap,DELTAx,Lpala,mac,Supaero,Uinf)
# dataout=[RHOelem,deltaCpelem,Aelem,VUelem,VLelem]
# elem (:,5:7) R punto de control
# elem (:,8:10) vector normal
# elem (:,11:13) V punto de control

# #Calculo par�metros aerodin�micos
#Cp=1-(V/Uinf)²-2/(Uinf²)*dfi/dt ver tesis SP p66
Pdin=0.5*rho*(VC)^2; # rho aire estandar (1.225 kg/m³ 0.0023769 slug/ft³). Es por VC por la adimensionalizacion, se cancela Uinf adimensional del DCp con la de Pdin y solo queda VC².

#salto de presi�n adimensional cp
cpelem=dataout(:,5); #Adimensional coeficiente de presión en punto de control
Aelem=dataout(:,6)*LC^2; #dimensional area del anillo

#normales en sistema B
nxelem=elem(:,8); #componente x de vector normal en punto de control
nyelem=elem(:,9); #componente y de vector normal en punto de control
nzelem=elem(:,10); #componente z de vector normal en punto de control

xcp=elem(:,5)*LC; #dimensional posicion x de punto de control
ycp=elem(:,6)*LC; #dimensional posicion y de punto de control
zcp=elem(:,7)*LC; #dimensional posicion z de punto de control

Rmom=cross(elem(:,5:7)*LC,elem(:,8:10),2); #M = rcp cross normal
Xmom=Rmom(:,1); #Momento torsi�n, dimensional
Ymom=Rmom(:,2); #Momento en eje rotor, dimensional

#######################################################
#Comprobaci�n aerodin�mica
#Amalla=sum(Aelem)
Aref=Supaero;

%calculo dx para primer fila de elementos
nodoizq=0.5*(nodosB(elem(1:n,1),1)+nodosB(elem(1:n,4),1))*LC; #promedio de nodos i j
nododer=0.5*(nodosB(elem(1:n,2),1)+nodosB(elem(1:n,3),1))*LC; #promedio de nodos k l
dxelem=(nododer-nodoizq)';

dzelem=Aelem./dxelem;
for i=1:n #barrido por perfil
  elementos=[i:n:(m-1)*n+i]; #identifico los elementos de la cuerda local.
%  nodoBA=0.5*(nodosB(elem(elementos,1),3)+nodosB(elem(elementos,2),3))*LC; #promedio de nodos i j
%  nodoBF=0.5*(nodosB(elem(elementos,3),3)+nodosB(elem(elementos,4),3))*LC; #promedio de nodos k l
%  dzelem=nodoBF-nodoBA;
  dzelem=Aelem(elementos)/dxelem(i); %dx es el mismo para toda la estacion
  FY2D(i)=Pdin*sum(nyelem(elementos).*cpelem(elementos).*dzelem);#fuerza normal de la seccion / dx
  FZ2D(i)=Pdin*sum(nzelem(elementos).*cpelem(elementos).*dzelem);#fuerza tangencial de la seccion / dx
  MX2D(i)=Pdin*sum(Xmom(elementos).*cpelem(elementos).*dzelem);#Momento en eje de pitch de la seccion / dx
%  FY2D(i)=Pdin*sum(nyelem(elementos).*cpelem(elementos).*dzelem);#fuerza normal de la seccion / dx
%  FZ2D(i)=Pdin*sum(nzelem(elementos).*cpelem(elementos).*dzelem);#fuerza tangencial de la seccion / dx
%  MX2D(i)=Pdin*sum(Xmom(elementos).*cpelem(elementos).*dzelem);#Momento en eje de pitch de la seccion / dx
endfor
%nodoizq=0.5*(nodosB(elem(1:n,1),1)+nodosB(elem(1:n,4),1))*LC; #promedio de nodos i j
%nododer=0.5*(nodosB(elem(1:n,2),1)+nodosB(elem(1:n,3),1))*LC; #promedio de nodos k l
%dxelem=(nododer-nodoizq)';

FY3D=sum(FY2D.*dxelem); #Fuerza total normal al plano de rotacion
FZ3D=sum(FZ2D.*dxelem); #Fuerza total dentro del plano de rotacion
MX3D=sum(MX2D.*dxelem);

%anillosBF=(m-1)*n+1:m*n;
%CMhinge=cpelem(anillosBF)/(2*m^2);

#######################################################
#########Cargas generalizadas##########################

xaero2=nodosB(elem(:,2),1)*LC;#xnodo j
xaero1=nodosB(elem(:,1),1)*LC;#xnodo i
deltax=xaero2-xaero1;
deltaz2=(nodosB(elem(:,3),3)-nodosB(elem(:,2),3))*LC; #znodok - znodoj
deltaz1=(nodosB(elem(:,4),3)-nodosB(elem(:,1),3))*LC; #znodol - znodoi
zmedio2=(nodosB(elem(:,3),3)+nodosB(elem(:,2),3))/2*LC; #promedio
zmedio1=(nodosB(elem(:,4),3)+nodosB(elem(:,1),3))/2*LC; #promedio
D=nxelem.*xcp+nyelem.*ycp+nzelem.*zcp;
#2 Puntos de gauss por cada elemento dimensionales!

%xg1=(xaero1+(1-1/sqrt(3))/2*deltax);
%xg2=(xaero1+(1+1/sqrt(3))/2*deltax);

xg1=(xaero1+(1-sqrt(3/5))/2*deltax);
xg2=(xaero1+(1+0)/2*deltax);
xg3=(xaero1+(1+sqrt(3/5))/2*deltax);

#Evaluacion de integrando en puntos de gauss
[uig1,uipg1,uippg1,vig1,vipg1,vippg1,wig1,wipg1,wippg1,tig1,tipg1,tippg1,nax,nflap,nedge,ntor]=modoseval_r1(xg1,L,nax,nflap);
[uig2,uipg2,uippg2,vig2,vipg2,vippg2,wig2,wipg2,wippg2,tig2,tipg2,tippg2,nax,nflap,nedge,ntor]=modoseval_r1(xg2,L,nax,nflap);
[uig3,uipg3,uippg3,vig3,vipg3,vippg3,wig3,wipg3,wippg3,tig3,tipg3,tippg3,nax,nflap,nedge,ntor]=modoseval_r1(xg3,L,nax,nflap);


pendiente=(deltaz2-deltaz1)./(xaero2-xaero1);
deltazg1=(deltaz1+pendiente.*(xg1-xaero1))';
deltazg2=(deltaz1+pendiente.*(xg2-xaero1))';
deltazg3=(deltaz1+pendiente.*(xg3-xaero1))';

pendiente=(zmedio2-zmedio1)./(xaero2-xaero1);
zmediog1=(zmedio1+pendiente.*(xg1-xaero1))';
zmediog2=(zmedio1+pendiente.*(xg2-xaero1))';
zmediog3=(zmedio1+pendiente.*(xg3-xaero1))';

xg1=xg1';
xg2=xg2';
xg3=xg3';

D=D';
nxelem=nxelem';
nyelem=nyelem';
nzelem=nzelem';
deltax=deltax';
cpelem=cpelem';

for i=1:nax
  Iu(i,:)=(5/9*deltazg1.*uig1(i,:)+8/9*deltazg2.*uig2(i,:)+5/9*deltazg3.*uig3(i,:)).*(deltax/2)./nyelem;

  Itg1(1,:)=(nzelem./(nyelem.^2).*(D-nxelem.*xg1)-((nzelem.^2)./(nyelem.^2)+1).*zmediog1.*deltazg1.*tig1(i,:));
  Itg2(1,:)=(nzelem./(nyelem.^2).*(D-nxelem.*xg2)-((nzelem.^2)./(nyelem.^2)+1).*zmediog2.*deltazg2.*tig2(i,:));
  Itg3(1,:)=(nzelem./(nyelem.^2).*(D-nxelem.*xg3)-((nzelem.^2)./(nyelem.^2)+1).*zmediog3.*deltazg3.*tig3(i,:));
  It(i,:)=(5/9*Itg1+8/9*Itg2+5/9*Itg3).*(deltax/2);
endfor

for i=1:nflap
  Iv1(i,:)=(5/9*deltazg1.*vig1(i,:)+8/9*deltazg2.*vig2(i,:)+5/9*deltazg3.*vig3(i,:)).*(deltax/2)./nyelem;
  Iw1(i,:)=(5/9*deltazg1.*wig1(i,:)+8/9*deltazg2.*wig2(i,:)+5/9*deltazg3.*wig3(i,:)).*(deltax/2)./nyelem;

  Iv2g1=deltazg1.*vipg1(i,:).*(D-nxelem.*xg1-nzelem.*zmediog1);
  Iv2g2=deltazg2.*vipg2(i,:).*(D-nxelem.*xg2-nzelem.*zmediog2);
  Iv2g3=deltazg3.*vipg3(i,:).*(D-nxelem.*xg3-nzelem.*zmediog3);
  Iv2(i,:)=(5/9*Iv2g1+8/9*Iv2g2+5/9*Iv2g3).*(deltax/2)./(nyelem.^2);
  
  Iw2g1=deltazg1.*zmediog1.*wipg1(i,:);
  Iw2g2=deltazg2.*zmediog2.*wipg2(i,:);
  Iw2g3=deltazg3.*zmediog3.*wipg3(i,:);
  Iw2(i,:)=(5/9*Iw2g1+8/9*Iw2g2+5/9*Iw2g3).*(deltax/2)./nyelem;
 
endfor

#Suma de todos los elementos de la superficie
Qu=sum(cpelem.*nxelem.*Iu,2);
Qv=sum(cpelem.*(nyelem.*Iv1-nxelem.*Iv2),2);
Qw=sum(cpelem.*(nzelem.*Iw1-nxelem.*Iw2),2); 
Qt=sum(cpelem.*It,2);

#Ensamblo vector de cargas generalizadas Qv,Qw,Qt
%Pdin=0.5*rho*(Uinf)^2; %Es Uinf porque el cp en Bernoulli se define respecto a la corriente de aire libre
%Faero=Pdin*[Qv';Qt']; #Modelo flexion-torsion Fung
Faero=Pdin*[Qv;Qw;Qt]; #Modelo flexion-torsion Pala 3D
%Faero=Pdin*[Qv;Qw;Qt]; #Modelo flexion-flexion-torsion Pala 3D

endfunction
