function [Ktg,Kg]=Kgeometrica3(x1,x2,qk,uip1,uip2,vip1,vip2,wip1,wip2,vipp1,vipp2,wipp1,wipp2,tip1,tip2,EA1,EA2,deltax,Ke,Kc,wy,nax,nflap,ntor)
## -*- texinfo -*-
## Rev1: 8/8/2014 modificado según revisión de las ecuaciones en "Modelo estructural nolineal de vigas tridimensionales indexado con energia potencial revisado 2.odt"
## Rev2: 10/9/2014 obtenido por solucion incremental iterativa
## Rev3: 9/8/2018 nax,nflap,ntor son parametros de entrada

## Author: NGT

%nax=size(uip1,1);
%nflap=size(vip1,1);
%ntor=size(tip1,1);

#Evaluo los terminos nolineales
	
	#Inicializo variables
	dudxk1=zeros(1,size(uip1,2));
	dudxk2=dudxk1;
	dvdxk1=zeros(1,size(vip1,2)); #thetaz
	dvdxk2=dvdxk1;
	dwdxk1=zeros(1,size(wip1,2));	#-thetay
	dwdxk2=dwdxk1;
	Ktg=zeros(nax+2*nflap+ntor,nax+2*nflap+ntor);
	Kg=zeros(nax+2*nflap+ntor,nax+2*nflap+ntor);
	
	#Calculo u0'
	for i=1:nax
		dudxk1=dudxk1+qk(i,1)*uip1(i,:);	
		dudxk2=dudxk2+qk(i,1)*uip2(i,:);	
  	endfor
    
	#Calculo v0'y w0'
	for i=1:nflap
		dvdxk1=dvdxk1+qk(nax+i,1)*vip1(i,:);
		dvdxk2=dvdxk2+qk(nax+i,1)*vip2(i,:);		
		dwdxk1=dwdxk1+qk(nax+nflap+i,1)*wip1(i,:);
		dwdxk2=dwdxk2+qk(nax+nflap+i,1)*wip2(i,:);		  
	endfor	
	
	
#Evaluo los coeficientes de la matriz Ktg
	
	#Terminos uv y uw
	#integro por gauss de 2 puntos: I = DX/2*(F1+F2)
	
	#Axial-flexion
	for i=1:nax
		for j=1:nflap
			Ktg(i,nax+j)=deltax/2*sum((EA1.*dvdxk1.*uip1(i,:).*vip1(j,:)+EA2.*dvdxk2.*uip2(i,:).*vip2(j,:)));#Kguv
			Ktg(i,nax+nflap+j)=deltax/2*sum((EA1.*dwdxk1.*uip1(i,:).*wip1(j,:)+EA2.*dwdxk2.*uip2(i,:).*wip2(j,:))); #Kguw
		endfor
	endfor
	Ktg(nax+1:nax+2*nflap,1:nax)=Ktg(1:nax,nax+1:nax+2*nflap)';#Kgvu,Kgwu
# 
	#Flexion-flexion
	for i=1:nflap
		for j=1:nflap
			Ktg(nax+i,nax+j)=deltax/2*sum((EA1.*(dudxk1+3/2*(dvdxk1.^2)+1/2*(dwdxk1.^2)).*vip1(i,:).*vip1(j,:)+EA2.*(dudxk2+3/2*(dvdxk2.^2)+1/2*(dwdxk2.^2)).*vip2(i,:).*vip2(j,:))); #KGvv
			Ktg(nax+i,nax+nflap+j)=deltax/2*sum((EA1.*dvdxk1.*dwdxk1.*vip1(i,:).*wip1(j,:)+EA2.*dvdxk2.*dwdxk2.*vip2(i,:).*wip2(j,:))); #KGvw
			Ktg(nax+nflap+i,nax+nflap+j)=deltax/2*sum((EA1.*(dudxk1+1/2*(dvdxk1.^2)+3/2*(dwdxk1.^2)).*wip1(i,:).*wip1(j,:)+EA2.*(dudxk2+1/2*(dvdxk2.^2)+3/2*(dwdxk2.^2)).*wip2(i,:).*wip2(j,:)));#KGww
		endfor
	endfor
	Ktg(nax+nflap+1:nax+2*nflap,nax+1:nflap)=Ktg(nax+1:nflap,nax+nflap+1:nax+2*nflap)';
###########################################################################################################################
#Evaluo el vector de fuerzas internas (para esquema incremental iterativo)

	#Axial-flexion

	for i=1:nax
		for j=1:nflap
			Kg(i,nax+j)=deltax/2*sum(((EA1/2).*dvdxk1.*uip1(i,:).*vip1(j,:)+(EA2/2).*dvdxk2.*uip2(i,:).*vip2(j,:)));#Kguv
			Kg(i,nax+nflap+j)=deltax/2*sum(((EA1/2).*dwdxk1.*uip1(i,:).*wip1(j,:)+(EA2/2).*dwdxk2.*uip2(i,:).*wip2(j,:))); #Kguw
		endfor
	endfor
	Kg(nax+1:nax+2*nflap,1:nax)=Kg(1:nax,nax+1:nax+2*nflap)';#Kgvu,Kgwu
# 
	#Flexion-flexion
	for i=1:nflap
		for j=1:nflap
			Kg(nax+i,nax+j)=deltax/2*sum(((EA1/2).*(dudxk1+dvdxk1.^2+dwdxk1.^2).*vip1(i,:).*vip1(j,:)+(EA2/2).*(dudxk2+dvdxk2.^2+dwdxk2.^2).*vip2(i,:).*vip2(j,:))); #KGvv
			Kg(nax+nflap+i,nax+nflap+j)=deltax/2*sum(((EA1/2).*(dudxk1+dvdxk1.^2+dwdxk1.^2).*wip1(i,:).*wip1(j,:)+(EA2/2).*(dudxk2+dvdxk2.^2+dwdxk2.^2).*wip2(i,:).*wip2(j,:))); #KGvv
		endfor
	endfor
	

endfunction

