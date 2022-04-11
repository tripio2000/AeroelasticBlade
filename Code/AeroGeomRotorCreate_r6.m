function [nodos,elem]=AeroGeomRotorCreate_r6(pala)
## -*- texinfo -*-
## @deftypefn {Function File} {} geomrotor (m,n,Ynodo,cuerda,Z,LC,VC,paso,elemflap,qflap)
## Entradas:
## m,n: divisiones de la superficie en cuerda y envergadura respectivamente
## Ynodo: vector de posiciones axiales 
## cuerda: cuerda en cada posicion Ynodo
## Z: vector de estado [q;qdot]
## LC: longitud de referencia mac/m
## VC: velocidad de referencia [m/s]
## paso: pitch de la pala en [rad]
## elemflap: Defino los elementos afectados por el flap
## qflap: coordenada generalizada asociada a la flexion del flap
## @end deftypefn
## Author: NGT
## Cuidado, el twist de los datos esta definido positivo antihorario
## Rev0: 01/09/2014 - Generaci�n de malla pala flexible en sistema rotante.
## Rev1: 31/08/2015 - Corro mallado 0.25 c/m hacia BFuga
## Rev2: 14/09/2015 - Uso Zpitch para ubicar perfiles. Zpitch se mide desde BAtaque hasta centro de pitch y se expresa en fracci�n de cuerda local.
## Rev3: 16/09/2015 - Incluyo comba de perfiles en variable airfoilnumber[afnumber] airfoildata[x,ymedia_1,ymedia_2,...,ymedia_n]
## Rev4: 16/12/2018 - Retomo posicion en c/4, mejoro comentarios
## Rev5: 18/04/2019 - Incluyo calculo de puntos de control en elem
## Rev6: 14/06/2019 - Uso estructura de datos para pala

m=pala.m;
n=pala.n;
raero=pala.raero;
cuerda=pala.cuerda;
twist=pala.twist;
zpitch=pala.zpitch;
airfoilnumber=pala.airfoilnumber;
airfoildata=pala.airfoildata;

#Creo anillos en la superficie
##################################################
	#armo tabla de conectividad
		k=1;
		for i=1:m
			for j=1:n
				ni=(i-1)*(n+1)+j;
				nj=(i-1)*(n+1)+j+1;
				nk=(i)*(n+1)+j+1;
				nl=(i)*(n+1)+j;
				elem(k,1:4)=[ni,nj,nk,nl];
				elem(k,14:15)=[i,j]; #i:z j:x
				k=k+1;
			endfor
		endfor

	
# Creo nodos de anillos adheridos siguiendo vista en planta
# los anillos están desfazados respecto de los paneles en cpanel/4
# m es la cantidad de paneles en cuerda y n en envergadura
##################################################
# Barro superficie zigzag axial - cuerda. i es en la cuerda, j es en la envergadura
	k=1;
	for i=1:m+1 %barro en direccion de la cuerda, m anillos tienen m+1 segmentos
		for j=1:n+1 %barro en envergadura, n anillos tienen n+1 segmentos
			zLE=-zpitch(j,1); %zpitch es la posicion del borde de ataque en fraccion de cuerda
%      zj(k,1)=(zLE+(i-1)/m)*cuerda(j,1);#;+0.25*cuerda(j,1)/m; 
      zj(k,1)=(zLE+(0.25+(i-1))/m)*cuerda(j,1);#coloco vortice en cpanel/4 para cumplir cond kutta
			airfoilcol=airfoilnumber(j,1); #veo que perfil es
      xperfil=airfoildata(:,1);
      yperfil=airfoildata(:,airfoilcol);
			yj(k,1)=interp1(xperfil,yperfil,(0.25+(i-1))/m,'extrap')*cuerda(j,1);
			k++;     
		endfor
	endfor
  

	#xj,yj,h,hdot,theta,thetadot son ((n+1)*(m+1),1)

	#Transformo vectores (n+1,1) en ((n+1)*(m+1))
	xj=repmat(raero,m+1,1);
	twist2=repmat(twist,m+1,1);
	mt=cos(-twist2);#el twist de los datos esta definido positivo antihorario
	nt=sin(-twist2);#el twist de los datos esta definido positivo antihorario

# Actualizo nodos por alabeo aerodinamico (twist)
##################################################
	#R0=[xj,yj,zj]; #sin alabeo
	#R0=[xj,zj.*nt,zj.*mt]; #con alabeo sin comba
	R0=[xj,-zj.*nt+yj.*mt,zj.*mt+yj.*nt];
	
	nodos(:,1:3)=R0;
  
  
     # Posicion del punto de control
	Ri=nodos(elem(:,1),1:3);
	Rj=nodos(elem(:,2),1:3);
	Rk=nodos(elem(:,3),1:3);
	Rl=nodos(elem(:,4),1:3);
	elem(:,5:7)=0.5*(Ri+0.5*(Rk-Ri)+Rj+0.5*(Rl-Rj)); 
	
endfunction
