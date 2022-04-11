function [Faero,FY3D,FY2D,FZ3D,FZ2D,MX3D,MX2D,nodosB,elem,nodoswakeB,G,dataout,sistema,placa]=AeroEvalFun_r1(dataout,pala,R0,elem,Z,d2qdt2,nax,nflap,ntor,wy,LC,TC,desprendimiento,nodoswakeB,conwake,VlibreB,arranque,G_old,deltaT,deltaTold,nodosB,VC,Uinf,rhoaire)
## Rev1: 14/6/19 uso estructura de datos pala
n=pala.n;
m=pala.m;
L=pala.L;
h=pala.h;
numelm=pala.numelm;
mac=pala.mac;
Supaero=pala.Supaero;

  
  #Actualización de geometría sábana adherida por desplazamiento y velocidad estructural
  [nodosB,elem]=AeroGeomRotorUpdate_r6(n,m,R0,elem,[Z(1:2*nflap+nax);Z(2*nflap+nax+1:end)/TC],nax,nflap,ntor,wy,LC,TC,L);#Z debe ser dimensional

  #Corrección de nodos de estela adheridos 
  #Determino nodos de borde libres
    if desprendimiento==1
		  nodosbordefuga=[m*(n+1)+1:(m+1)*(n+1)];#con esquinas
      nodosstartvrtx=[nodosbordefuga];
	  elseif desprendimiento==2
      nodosbordeder=[n+1:n+1:(m+1)*(n+1)];
      nodosbordefuga=[m*(n+1)+1:(m+1)*(n+1)-1];#sin esquina derecha
      nodosstartvrtx=[nodosbordeder,nodosbordefuga];
    elseif desprendimiento==3
      nodosbordeizq=[1:n+1:m*(n+1)+1];
  		nodosbordeder=[n+1:n+1:(m+1)*(n+1)];
  		nodosbordefuga=[m*(n+1)+2:(m+1)*(n+1)-1];#sin esquinas
  		nodosstartvrtx=[nodosbordeizq,nodosbordeder,nodosbordefuga];
    endif
    #Determino posiciones adheridas
    nnodoslibres=length(nodosstartvrtx);
    #corrijo posicion de nodos adheridos
    nodoswakeB(end-nnodoslibres+1:end,1:3)=nodosB(nodosstartvrtx,1:3);
  
  #Evaluación de circulaciones 
  Ri=nodoswakeB(conwake(:,1),1:3);
	Rj=nodoswakeB(conwake(:,2),1:3);
	Rk=nodoswakeB(conwake(:,3),1:3);
	Rl=nodoswakeB(conwake(:,4),1:3);
	Gwake=conwake(:,5);
	[G,Vindvecfull]=AeroSolPHI_rev4(nodosB,elem,VlibreB,Ri,Rj,Rk,Rl,Gwake,numelm,arranque); %calculo circulaciones en t+dt
  DGDt=(G-G_old)/deltaTold;

  #Evaluación de salto de presión
	[dataout]=AeroPressureSP_rev11(nodosB,elem,numelm,VlibreB,G,DGDt,nodoswakeB,conwake,m,n,Vindvecfull,desprendimiento); 

  #Evaluación de cargas
  [Faero,FY3D,FY2D,FZ3D,FZ2D,MX3D,MX2D]=AeroLoads_rev10(VC,LC,dataout,elem,nodosB,mac,Supaero,Uinf,L,rhoaire,wy,nax,nflap,n,m);

  endfunction
  