function [v0,w0,thetax,thetay,thetaz,dotv0,dotw0,dotthetax,dotthetay,dotthetaz]=desplazamientos_r4(xb,L,Z,nax,nflap)
## Calculo de desplazamientos del centro elástico para viga con funciones de forma modales
## Z es el vector de estado [qi;qdoti] Zdot es el vector de estado [qdoti;qddoti]
## Rev1: 3/9/15 incluyo velocidades
## Rev2: 15/10/15 elimino grado axial
## Rev3: 22/12/15 agrego aceleraciones
## Rev4: 9/2/19 elimino aceleraciones
	#Evaluacion de ff forma en puntos xb. xb es vector fila xb[1,:]
  [ui,uip,uipp,vi,vip,vipp,wi,wip,wipp,ti,tip,tipp,nax,nflap,nedge,ntor]=modoseval_r1(xb,L,nax,nflap);
  #Las salidas son [ui,x] donde i=1:nmodos

	#Desplazamientos
	qv=Z(1:nflap,:);
	qw=Z(nflap+1:2*nflap,:);
	qtx=Z(2*nflap+1:2*nflap+nax,:);	

	v0=vi'*qv; #vi * qv
	w0=wi'*qw;#wi * qw
	thetax=ti'*qtx;#ti * qt
	thetay=vip'*qv;#ti * qt
	thetaz=wip'*qw;#ti * qt

	qend=2*nflap+nax;

	#Velocidades
	dotqv=Z(qend+1:qend+nflap,:);
	dotqw=Z(qend+nflap+1:qend+2*nflap,:);
	dotqtx=Z(qend+2*nflap+1:end,:);
	
	dotv0=vi'*dotqv; #vi * dot qv
	dotw0=wi'*dotqw; #wi * dot qw
	dotthetax=ti'*dotqtx; #ti * dot qt
	dotthetay=vip'*dotqv; #thetay=dv0/dx, vip dot qv
	dotthetaz=wip'*dotqw; #thetaz=dw0/dx, #wip dot qw 

%	#Aceleraciones
%	ddotqv=Zdot(qend+1:qend+nflap,:);
%	ddotqw=Zdot(qend+nflap+1:qend+2*nflap,:);
%	ddotqtx=Zdot(qend+2*nflap+1:end,:);
%	
%	ddotv0=vi'*ddotqv;
%	ddotw0=wi'*ddotqw;
%	ddotthetax=ti'*ddotqtx;

endfunction


