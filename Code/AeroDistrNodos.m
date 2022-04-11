function [raero,n]=AeroDistrNodos(espaciado,raeroimport,n,m)
#Metodo de posicionado de nodos en envergadura:0 uniforme, 1 refinamiento senoidal hacia la punta, 2 refinamiento senoidal, 2 refinamiento senoidal incompleto.
# n es el numero de divisiones, tengo n+1 estaciones
if espaciado ==0
	#Espaciado uniforme
	raero=linspace(raeroimport(1),raeroimport(end),n+1)';
elseif espaciado ==1
	#Refinamiento punta senoidal
  thetatip=0.4*pi;
	dthetaaero=(thetatip)/(n);
	thetaaero=0:dthetaaero:thetatip;
	x0=raeroimport(1);
	xtip=raeroimport(end);
	raero=x0+(xtip-x0)*sin(thetaaero)/sin(thetatip);
	raero=raero';
elseif espaciado==2
	#Refinamiento punta y raiz senoidal 
  thetahub=-0.5*pi;
  thetatip=0.5*pi;
	dthetaaero=(thetatip-thetahub)/n;
	thetaaero=thetahub:dthetaaero:thetatip;
	Laero=raeroimport(end)-raeroimport(1);
	x0=raeroimport(1);
	xm=x0+Laero/2;
	raero=xm+Laero/2*sin(thetaaero)/sin(thetatip);
	raero=raero';
elseif espaciado==4
	#cuadrados pala WindPACT
  r=raeroimport(1);
  i=1;
  while r<raeroimport(end)-0.75*0.9/m
    raero(i)=r;
    if r<0.25*35
      cuerdalocal=2.389+(2.8-2.389)/(8.75-4.374)*(r-raeroimport(1));
    else
      cuerdalocal=2.8+(0.906-2.8)/(35-8.75)*(r-8.75);
    endif
    dz=cuerdalocal/m;
    dr=dz;
    r=r+dr;
    i++;
  endwhile
  raero(i)=raeroimport(end);
  raero=raero';
  n=i-1;
   
endif

