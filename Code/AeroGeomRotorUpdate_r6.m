function [nodos,elem]=AeroGeomRotorUpdate_r6(n,m,R0,elem,Z,nax,nflap,ntor,wy,LC,TC,L)
## Rev1: Uso funciones de forma modales
## Rev2: Incluyo deformacion flap, elimino desplazamiento elastico axial
## Rev3: Actualizo desplazamiento por flap
## Rev5: Uso estructura placa
## Rev6: quito flap

## Z:vector de coordenadas generalizadas
## El Z de entrada es adimensional por el tiempo, R0 es dimensional
## R0: posiciones de los nodos de la pala sin deformar 
## Los nodos aerodin�micos se barren en zig-zag siguiendo direcci�n axial y luego cuerda. Las posiciones axiales se obtienen de xaero. La pala tiene m+1 estaciones axiales y n+1 puntos en la cuerda.
## Author: NGT

##CUIDADO la geometria se actualiza en dimensional!!!!

################################################################
#Desplazamientos del centro el�stico de todos los nodos de la grilla (R0) medidos desde el hub 
%	[v0,w0,thetax,thetay,thetaz,dotv0,dotw0,dotthetax,dotthetay,dotthetaz]=desplazamientos_fung_r4((R0(:,1))',L,Z);
  [v0,w0,thetax,thetay,thetaz,dotv0,dotw0,dotthetax,dotthetay,dotthetaz]=desplazamientos_r4((R0(:,1))',L,Z,nax,nflap); #Los resultados son vectores columna. v0[:,1]
  
  #Desplazamientos de nodos de grilla aerodin�mica
	Relastico=[(R0(:,3).*thetay)-(R0(:,2).*thetaz), v0-(R0(:,3).*thetax), w0+(R0(:,2).*thetax)];
	#Relastico=[u0'+(R0(:,3).*thetay')-(R0(:,2).*thetaz'), v0'-(R0(:,3).*thetax'), w0'+(R0(:,2).*thetax')];
  
	#Actualizo posicion por deformacion el�stica:
	nodos(:,1:3)=R0+Relastico;
  
################################################################
# Velocidades

	#Velocidad por deformacion el�stica 	#vel=[udot;vdot;wdot]+[0,-wz,wy;wz,0,-wx;-wy,wx,0]*(ro+u)
	Velastico=[(R0(:,3).*dotthetay)-(R0(:,2).*dotthetaz), dotv0-(R0(:,3).*dotthetax), dotw0+(R0(:,2).*dotthetax)];
  
	#Actualizo velocidad por deformacion el�stica:
	nodos(:,4:6)=Velastico;

	#Actualizo velocidades por marco rotante
	#Velocidad absoluta en sistema B V = udot + wy*[0,0,1;0,0,0;-1,0,0] cdot (R0+u) # wy*[0,0,1;0,0,0;-1,0,0] cdot (R0+u) = [wy (R0+u)_z;0;-wy(R0+u)_x]
	nodos(:,4:6)=nodos(:,4:6)+wy*[nodos(:,3),zeros((n+1)*(m+1),1),-nodos(:,1)];

	
################################################################
# Adimensionalizacion
	#Adimensionalizo posiciones de los nodos:
	nodos(:,1:3)=nodos(:,1:3)/LC;
	#Adimensionalizo velocidades de los nodos:
	nodos(:,4:6)=nodos(:,4:6)*(TC/LC);

################################################################
# Determino posicion y velocidades de los elementos

  # Posicion del punto de control
	Ri=nodos(elem(:,1),1:3);
	Rj=nodos(elem(:,2),1:3);
	Rk=nodos(elem(:,3),1:3);
	Rl=nodos(elem(:,4),1:3);
	elem(:,5:7)=0.5*(Ri+0.5*(Rk-Ri)+Rj+0.5*(Rl-Rj));
	
	# Vector normal 
	normales=cross(Ri-Rk,Rl-Rj);
	divisor=sqrt(normales(:,1).^2+normales(:,2).^2+normales(:,3).^2);
	elem(:,8:10)=normales./divisor;
	
	# Velocidad del punto de control
	Vi=nodos(elem(:,1),4:6);
	Vj=nodos(elem(:,2),4:6);
	Vk=nodos(elem(:,3),4:6);
	Vl=nodos(elem(:,4),4:6);
	elem(:,11:13)=0.25*(Vi+Vj+Vk+Vl); #Promedio de las velocidades nodales (adimensionales)

endfunction
