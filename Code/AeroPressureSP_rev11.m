function [dataout]=AeroPressureSP_rev11 (nodosB,elem,numelm,VlibreB,G,DGDt,nodoswakeB,conwake,m,n,Vindwakevec,desprendimiento)
## elimino anillos desprendidos por borde izquierdo segun convectar_rev10
## optimizo codigo
## Actualizo Vinducida_rev2
## Rev8: 18/8/2015 Corrijo Cp
## Rev10: 16/12/18 Agrego DCp por Kutta Joukowski. Fracaso.
## Rev11: 17/12/18 Agrego selecci�n de bordes libres para rhovec
	
#Determino tama�os de variables
numwake=rows(conwake); #cantidad de anillos en estela
%tamwake=(m+n); #anillos convectados por paso de tiempo si hay estela 1 lateral
tamwake=(n); #anillos convectados por paso de tiempo si hay estela solo fuga
nconvect=numwake/tamwake; #n�mero de estelas convectadas

#Asigno variables
ni=nodosB(elem(:,1),1:3);
nj=nodosB(elem(:,2),1:3);
nk=nodosB(elem(:,3),1:3);
nl=nodosB(elem(:,4),1:3);
wi=nodoswakeB(conwake(:,1),1:3);
wj=nodoswakeB(conwake(:,2),1:3);
wk=nodoswakeB(conwake(:,3),1:3);
wl=nodoswakeB(conwake(:,4),1:3);
Gw=conwake(:,5);

#Determino velocidades sobre nodos de superficie en el sistema B
#la velocidad inducida por la estela ya es conocida y es un arg de entrada
##########################################################################
#Velocidad inducida por superficie
[Vindsurfvec]=AeroVinducida_rev5(elem(:,5:7),G,ni,nj,nk,nl,1);
%disp('Vbound en cp para cargas')
#Velocidad libre
VlibreBvec=repmat(VlibreB,numelm,1); #Es adimensional por VC=LC/TC!
#Velocidad por cinem�tica de la superficie (contempla vibraciones el�sticas y marco rotante)
VdinamicaB=elem(:,11:13); #OJO que las velocidades en elem son adimensionales por VC=LC/TC!

 
######Asigno circulaciones para cada segmento de cada anillo
#rho1: Oeste
#rho2: Norte 
#rho3: Sur
#rho4: Este
####################################################
%desprendimiento=1; #1: borde de fuga, 2: borde de fuga y puntera derecha, 3: borde de fuga y ambas punteras

if desprendimiento==1 
  wakeelem=n*(nconvect-1)+1:n*(nconvect-1)+n; #lista de elementos de estela si solo borde de fuga
  %borde ataque izquierdo (i==1)
  %'elemento:',i,'elemento:BAI'
  ibai=1;
  rho1vec(ibai,1)=G(ibai); 
  rho2vec(ibai,1)=2*G(ibai); #se multiplica por 2 porque es una superficie sustentadora y no se conecta con otro panel
  rho3vec(ibai,1)=-G(ibai)+G(ibai+n);
  rho4vec(ibai,1)=-G(ibai)+G(ibai+1);
  %borde ataque derecho (i==n)
  %'elemento:',i,'elemento:BAD';
  ibad=n;
  rho1vec(ibad,1)=G(ibad)-G(ibad-1); 
  rho2vec(ibad,1)=2*G(ibad); #se multiplica por 2 porque es una superficie sustentadora y no se conecta con otro panel
  rho3vec(ibad,1)=-G(ibad)+G(ibad+n);
  rho4vec(ibad,1)=-G(ibad);	
  %borde fuga izquierdo (i==(m-1)*n+1)
  %   'elemento:',i,'elemento:BFI';
  ibfi=(m-1)*n+1;
  rho1vec(ibfi,1)=G(ibfi);
  rho2vec(ibfi,1)=G(ibfi)-G(ibfi-n);
  rho3vec(ibfi,1)=-G(ibfi)+Gw(wakeelem(1),1); %si los nodos esquina estan en borde de fuga
  rho4vec(ibfi,1)=-G(ibfi)+G(ibfi+1);
  %borde fuga derecho (i==m*n)
  %  'elemento:',i,'elemento:BFD';
  ibfd=m*n;
  rho1vec(ibfd,1)=G(ibfd)-G(ibfd-1); 
  rho2vec(ibfd,1)=G(ibfd)-G(ibfd-n);
  rho3vec(ibfd,1)=-G(ibfd)+Gw(wakeelem(n),1); %si los nodos esquina estan en borde de fuga
  rho4vec(ibfd,1)=-G(ibfd);	%si no hay estela lateral
  %borde fuga (i>(m-1)*n+2 & i<m*n-1)
  %'elemento:',i,'elemento:BF';
  ibf=(m-1)*n+2:m*n-1;
  columna=2:n-1;
  rho1vec(ibf,1)=G(ibf)-G(ibf-1); 
  rho2vec(ibf,1)=G(ibf)-G(ibf-n);
  rho3vec(ibf,1)=-G(ibf)+Gw(wakeelem(columna),1); %si los nodos esquina estan en borde de fuga
  rho4vec(ibf,1)=-G(ibf)+G(ibf+1);
  %lateral izquierdo
  %'elemento:',i,'elemento:LI';
  ili=(n+1:n:n*(m-1));
  rho1vec(ili,1)=G(ili);
  rho2vec(ili,1)=G(ili)-G(ili-n);
  rho3vec(ili,1)=-G(ili)+G(ili+n);
  rho4vec(ili,1)=-G(ili)+G(ili+1);
  %lateral derecho
  %'elemento:',i,'elemento:LD';
  ild=(2*n:n:n*(m-1));
  rho1vec(ild,1)=G(ild)-G(ild-1); 
  rho2vec(ild,1)=G(ild)-G(ild-n);
  rho3vec(ild,1)=-G(ild)+G(ild+n);
  rho4vec(ild,1)=-G(ild);

  elseif desprendimiento==2 #borde de fuga y extremo derecho
 
  wakeelem=numwake-(n+m)+1:numwake; #lista de anillos de estela vecinos, ordenados por [m derechos, n-1 fuga, esquina]

  %borde ataque izquierdo (i==1)
  %'elemento:',i,'elemento:BAI'
  ibai=1;
  rho1vec(ibai,1)=G(ibai); 
  rho2vec(ibai,1)=2*G(ibai);#se multiplica por 2 porque es una superficie sustentadora y no se conecta con otro panel
  rho3vec(ibai,1)=-G(ibai)+G(ibai+n);
  rho4vec(ibai,1)=-G(ibai)+G(ibai+1);

  %borde ataque derecho (i==n)
  %'elemento:',i,'elemento:BAD';
  ibad=n;
  rho1vec(ibad,1)=G(ibad)-G(ibad-1); 
  rho2vec(ibad,1)=2*G(ibad);#se multiplica por 2 porque es una superficie sustentadora y no se conecta con otro panel
  rho3vec(ibad,1)=-G(ibad)+G(ibad+n);
  rho4vec(ibad,1)=-G(ibad)+Gw(wakeelem(1),1);#borde lateral
  
  %borde fuga izquierdo (i==(m-1)*n+1)
  %   'elemento:',i,'elemento:BFI';
  ibfi=(m-1)*n+1;
  rho1vec(ibfi,1)=G(ibfi);
  rho2vec(ibfi,1)=G(ibfi)-G(ibfi-n);
  rho3vec(ibfi,1)=-G(ibfi)+Gw(wakeelem(m+1),1);#borde fuga
  rho4vec(ibfi,1)=-G(ibfi)+G(ibfi+1);
	
  %borde fuga derecho (i==m*n)
  %  'elemento:',i,'elemento:BFD';
  ibfd=m*n;
  rho1vec(ibfd,1)=G(ibfd)-G(ibfd-1); 
  rho2vec(ibfd,1)=G(ibfd)-G(ibfd-n);
  rho3vec(ibfd,1)=-G(ibfd)+Gw(wakeelem(end),1); #borde fuga
  rho4vec(ibfd,1)=-G(ibfd)+Gw(wakeelem(m),1); #borde lateral
			
  %borde fuga (i>(m-1)*n+2 & i<m*n-1)
  %'elemento:',i,'elemento:BF';
  ibf=(m-1)*n+2:m*n-1;
  columna=m+1:m+n-2;
  rho1vec(ibf,1)=G(ibf)-G(ibf-1); 
  rho2vec(ibf,1)=G(ibf)-G(ibf-n);
  rho3vec(ibf,1)=-G(ibf)+Gw(wakeelem(columna),1); %si los nodos esquina estan en borde de fuga
  rho4vec(ibf,1)=-G(ibf)+G(ibf+1);

  %lateral izquierdo
  %'elemento:',i,'elemento:LI';
  ili=(n+1:n:n*(m-1));
  %fila=(ili-1)/n+1;
  rho1vec(ili,1)=G(ili);
  rho2vec(ili,1)=G(ili)-G(ili-n);
  rho3vec(ili,1)=-G(ili)+G(ili+n);
  rho4vec(ili,1)=-G(ili)+G(ili+1);
		
  %lateral derecho
  %'elemento:',i,'elemento:LD';
  ild=(2*n:n:n*(m-1));
  columna=2:m-1;
  rho1vec(ild,1)=G(ild)-G(ild-1); 
  rho2vec(ild,1)=G(ild)-G(ild-n);
  rho3vec(ild,1)=-G(ild)+G(ild+n);
  rho4vec(ild,1)=-G(ild)+Gw(wakeelem(columna),1);

  
  elseif desprendimiento==3
   wakeelem=(numwake-(n+2*m)+1):numwake; #lista de anillos de estela vecinos, ordenados por [m izquierdos, m derechos, n-2 fuga, 2 esquina]

  %borde ataque izquierdo (i==1)
  %'elemento:',i,'elemento:BAI'
  ibai=1;
  rho1vec(ibai,1)=G(ibai)-Gw(wakeelem(1),1); #contacto con anillo estela
  rho2vec(ibai,1)=2*G(ibai);#se multiplica por 2 porque es una superficie sustentadora y no se conecta con otro panel
  rho3vec(ibai,1)=-G(ibai)+G(ibai+n);
  rho4vec(ibai,1)=-G(ibai)+G(ibai+1);

  %borde ataque derecho (i==n)
  %'elemento:',i,'elemento:BAD';
  ibad=n;
  rho1vec(ibad,1)=G(ibad)-G(ibad-1); 
  rho2vec(ibad,1)=2*G(ibad);#se multiplica por 2 porque es una superficie sustentadora y no se conecta con otro panel
  rho3vec(ibad,1)=-G(ibad)+G(ibad+n);
  rho4vec(ibad,1)=-G(ibad)+Gw(wakeelem(m+1),1); #contacto con anillo estela
  
  %borde fuga izquierdo (i==(m-1)*n+1)
  %   'elemento:',i,'elemento:BFI';
  ibfi=(m-1)*n+1;
  rho1vec(ibfi,1)=G(ibfi)-Gw(wakeelem(m),1);
  rho2vec(ibfi,1)=G(ibfi)-G(ibfi-n);
  rho3vec(ibfi,1)=-G(ibfi)+Gw(wakeelem(end-1),1);#borde fuga
  rho4vec(ibfi,1)=-G(ibfi)+G(ibfi+1);
	
  %borde fuga derecho (i==m*n)
  %  'elemento:',i,'elemento:BFD';
  ibfd=m*n;
  rho1vec(ibfd,1)=G(ibfd)-G(ibfd-1); 
  rho2vec(ibfd,1)=G(ibfd)-G(ibfd-n);
  rho3vec(ibfd,1)=-G(ibfd)+Gw(wakeelem(end),1); #borde fuga
  rho4vec(ibfd,1)=-G(ibfd)+Gw(wakeelem(2*m),1); #borde lateral
			
  %borde fuga (i>(m-1)*n+2 & i<m*n-1)
  %'elemento:',i,'elemento:BF';
  ibf=(m-1)*n+2:m*n-1;
  columna=2*m+1:2*m+n-2;
  rho1vec(ibf,1)=G(ibf)-G(ibf-1); 
  rho2vec(ibf,1)=G(ibf)-G(ibf-n);
  rho3vec(ibf,1)=-G(ibf)+Gw(wakeelem(columna),1); %si los nodos esquina estan en borde de fuga
  rho4vec(ibf,1)=-G(ibf)+G(ibf+1);

  %lateral izquierdo
  %'elemento:',i,'elemento:LI';
  ili=(n+1:n:n*(m-1));
  columna=2:m-1;
  %fila=(ili-1)/n+1;
  rho1vec(ili,1)=G(ili)-Gw(wakeelem(columna),1);
  rho2vec(ili,1)=G(ili)-G(ili-n);
  rho3vec(ili,1)=-G(ili)+G(ili+n);
  rho4vec(ili,1)=-G(ili)+G(ili+1);
		
  %lateral derecho
  %'elemento:',i,'elemento:LD';
  ild=(2*n:n:n*(m-1));
  columna=m+2:2*m-1;
  rho1vec(ild,1)=G(ild)-G(ild-1); 
  rho2vec(ild,1)=G(ild)-G(ild-n);
  rho3vec(ild,1)=-G(ild)+G(ild+n);
  rho4vec(ild,1)=-G(ild)+Gw(wakeelem(columna),1);
  
  endif

%borde ataque (i>1 & i<n)
%	'elemento:',i,'elemento:BA';
iba=2:n-1;
rho1vec(iba,1)=G(iba)-G(iba-1); 
rho2vec(iba,1)=2*G(iba);#se multiplica por 2 porque es una superficie sustentadora y no se conecta con otro panel
rho3vec(iba,1)=-G(iba)+G(iba+n);
rho4vec(iba,1)=-G(iba)+G(iba+1);
  
% segmentos Interiores
iusados=union(union(union(ibai,ibad),union(iba,ibfi)),union(union(ibfd,ibf),union(ili,ild)));
%'elemento:',i,'elemento:INT';
iint=setdiff(1:n*m,iusados);
rho1vec(iint,1)=G(iint)-G(iint-1); 
rho2vec(iint,1)=G(iint)-G(iint-n);
rho3vec(iint,1)=-G(iint)+G(iint+n);
rho4vec(iint,1)=-G(iint)+G(iint+1);
rhovec=[rho1vec,rho2vec,rho3vec,rho4vec];

###################################################
		
#Asigno variables
rcpi=elem(:,5:7);
L1=ni-nl;
L2=nj-ni;
L3=nk-nl;
L4=nj-nk;
#Calculo areas
A1=cross(L1,L3,2);
A2=cross(L2,L4,2);
Ai(:,1)=0.5*sqrt(dot(A1,A1,2))+0.5*sqrt(dot(A2,A2,2));
#Armo vector de circulaciones del anillo
rhov=0.5*(rhovec(:,1).*L1+rhovec(:,2).*L2+rhovec(:,3).*L3+rhovec(:,4).*L4);
#Salto de velocidad tangencial por circulaci�n local
deltaV =-cross(elem(:,8:10),rhov(:,:) ,2 ) ./ Ai(:,1);
#Vm=VlibreBvec+Vindwakevec+VdinamicaB;
Vm=VlibreBvec+Vindwakevec;
#Proyecto en direccion tangente (Vtangente = Vtotal - Vnormal)
Vm=Vm-( dot(Vm,elem(:,8:10),2) .*elem(:,8:10) ); 
Vcp=VdinamicaB;
deltaCp= 2*( dot( Vm-Vcp ,deltaV, 2 )+ DGDt );

VU=Vm+deltaV/2;
VL=Vm-deltaV/2;

%#Salto de presion por teorema Kutta Joukowski
%deltax=norm(L2,2,'rows');
%rholocal=rho2vec;
%deltaCpKJ=2*rholocal./(deltax); %es el DCP*Qinf

dataout=[rhovec,deltaCp,Ai,VU,VL]; #rhovec son las circulaciones de cada segmento


endfunction

