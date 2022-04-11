function [Vindvec]=AeroVinducida_rev5(Rp,G,Ri,Rj,Rk,Rl,sumar)
## -*- texinfo -*-
## @deftypefn {Function File} {} Vinducida_rev ()
## rcp es el vector de los puntos de evaluaci�n
## ri,rj,rk,rl son vectores de nodos de los elementos emisores
## G es el vector de las vorticidades de cada emisor
## nrx=numelm; #tama�o receptor
## ntx=numwake #tama�o emisor
## Rev1 3/10/14 Limito tama�o asignable en memoria
## Rev2 8/9/15 uso voring_bsx
## Rev3 4/12/18 agrego logica para reducir calculos
## Rev4 4/12/18 Limito tama�o asignable en memoria
## Rev5 9/12/18 uso distancial radial al cp, en vez de d normal.

nrx=rows(Rp);
ntx=rows(G);

nrxmax=2000;#tamano del paquete de datos (1000 para 1.2Gb RAM)
if nrx<nrxmax
	nrxmax=nrx;
endif

ioxfer=ceil(nrx/nrxmax); #cantidad de transacciones de i/o, siempre tomo redondeo mayor
nrxi=zeros(1,nrxmax);

#Reordeno todos contra todos, el tamaño final es ntx*nrx

  
	for i=1:ioxfer
	# divido los nrx por chunks de tamaño nrxmax
  # Defino lista de nodos con un vector con indices nrxi
		if i==1
			nrxi=1:nrxmax; %chunk completo
      Gvec=repmat(G,nrxmax,1);
	    Rivec=repmat(Ri,nrxmax,1);
	    Rjvec=repmat(Rj,nrxmax,1);
	    Rkvec=repmat(Rk,nrxmax,1);
	    Rlvec=repmat(Rl,nrxmax,1);
		elseif i==ioxfer 	# #el ultimo es de tamaño menor
			pivot=nrxi(1,end)+1;
			nrxi=pivot:nrx; #desde donde quedo el chunk anterior hasta el final
			nrxmax=length(nrxi);
			Gvec=repmat(G,nrxmax,1);
			Rivec=repmat(Ri,nrxmax,1);
			Rjvec=repmat(Rj,nrxmax,1);
			Rkvec=repmat(Rk,nrxmax,1);
			Rlvec=repmat(Rl,nrxmax,1);	
		else
			nrxi=nrxi+nrxmax*ones(1,nrxmax); %chunk completo
      Gvec=repmat(G,nrxmax,1);
	    Rivec=repmat(Ri,nrxmax,1);
	    Rjvec=repmat(Rj,nrxmax,1);
	    Rkvec=repmat(Rk,nrxmax,1);
	    Rlvec=repmat(Rl,nrxmax,1);
		endif
    #Llevo a tama�o maxmemoria (repito columnas)
	  Rx=[repmat(Rp(nrxi,1),1,ntx)]';
	  Ry=[repmat(Rp(nrxi,2),1,ntx)]';
	  Rz=[repmat(Rp(nrxi,3),1,ntx)]';
	  #Traspongo para aplicar funci�n de forma vectorial
	  Rx=reshape(Rx,ntx*nrxmax,1);
	  Ry=reshape(Ry,ntx*nrxmax,1);
	  Rz=reshape(Rz,ntx*nrxmax,1);
    Vind=zeros(ntx*nrxmax,3);
 
 reducir_interacciones=0;
 
 if reducir_interacciones==1;   
    #Determino cuales interacciones no son despreciables
     Rcp=(Rivec+Rjvec+Rkvec+Rlvec)/4; #centro del anillo o control point (cp)
%     a=min([norm(Rivec,2,'rows'),norm(Rjvec,2,'rows'),norm(Rkvec,2,'rows'),norm(Rlvec,2,'rows')],[],2);
%     L1=Rkvec-Rivec;
%     L2=Rjvec-Rlvec;
%     a=(norm(L1,2,'rows')+norm(L2,2,'rows'))/2; #longitud del anillo
    L1=Rivec-Rlvec;
    L2=Rjvec-Rivec;
    L3=Rkvec-Rlvec;
    L4=Rjvec-Rkvec;
    #Calculo areas
    A1=cross(L1,L3,2);
    A2=cross(L2,L4,2);
    A=0.5*sqrt(dot(A1,A1,2))+0.5*sqrt(dot(A2,A2,2));
    a=sqrt(A);  #longitud aproximada del anillo asumiendo forma cuadrada
    nhat=cross(L1,L2,2);
    nhat=nhat./norm(nhat,2,'rows');
%     d=abs(dot([Rx,Ry,Rz]-Rcp,nhat,2)); %distancia normal del anillo al punto
     d=norm([Rx,Ry,Rz]-Rcp,2,'rows'); %distancia del punto a evaluar al cp del anillo
     anillos=find((d./a)<=20); #elijo en base a distancia/tamaño
     puntos=find((d./a)>20); #eligo en base a distancia/tamaño
    Vind(anillos,:)=AeroVorRing0_rev1(Gvec(anillos,:),Rivec(anillos,:),Rjvec(anillos,:),Rkvec(anillos,:),Rlvec(anillos,:),[Rx(anillos,:),Ry(anillos,:),Rz(anillos,:)]); #Calculo velocidad inducida
    Vind(puntos,:)=AeroDblPnt(Gvec(puntos,:),Rcp(puntos,:),nhat(puntos,:),A(puntos,:),[Rx(puntos,:),Ry(puntos,:),Rz(puntos,:)]);
  else
    Vind=AeroVorRing0_rev1(Gvec,Rivec,Rjvec,Rkvec,Rlvec,[Rx,Ry,Rz]); #Calculo velocidad inducida
  endif
  
    
    if sumar==0 #no combinar las velocidades inducidas individuales (para armar la matriz A)
     Vindvecx(:,nrxi)=[reshape(Vind(:,1),ntx,nrxmax)];
     Vindvecy(:,nrxi)=[reshape(Vind(:,2),ntx,nrxmax)];
     Vindvecz(:,nrxi)=[reshape(Vind(:,3),ntx,nrxmax)];
     Vindvec=[Vindvecx,Vindvecy,Vindvecz];
    else #combinar las velocidades inducidas individuales
     Vindvec(nrxi,:)=[sum(reshape(Vind(:,1),ntx,nrxmax),1)',sum(reshape(Vind(:,2),ntx,nrxmax),1)',sum(reshape(Vind(:,3),ntx,nrxmax),1)'];
    endif
    
endfor


endfunction
