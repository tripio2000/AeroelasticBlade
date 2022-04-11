function [nodoswake,conwake]=AeroConvect_rev14(nodos,elem,G,nodoswake,conwake,Vlibre,deltaT,m,n,initestela,wy,desprendimiento)
## Evaluo velocidades en estela y muevo por DT.
## Vectorizado conveccion de estela 30-4-2014
## Implemento l�gica para evaluar solamente campo cercano 5-5-2014
## Modularizo evaluaci�n de vinducida
## Reduzco uso de memoria
## Depuraci�n y ordenamiento
## Resuelvo todo en sistema B
## RK2 para movimiento por wy
## Rev9 elimino opciones de conveccion
## Rev10 no se desprende del borde Izquierdo (Raiz-hub)
## Rev11 agrego truncamiento de estela
## Rev12 actualizo vinducida_rev2
## Rev13 implemento aparte el esquema explicito para velocidad en estela
## Rev14 17/12/18 Agrego selecci�n de bordes libres

	#Los nodos, nodoswake y Vlibre estan en el sistema rotante (B)
  
    
	#Asigno variables
	numelm=size(G,1); 
	#Determino esquinas de cada elemento
	niB=nodos(elem(:,1),1:3); 
	njB=nodos(elem(:,2),1:3);
	nkB=nodos(elem(:,3),1:3);
	nlB=nodos(elem(:,4),1:3);	

	if initestela==0 #Ya hay estela, 
    rowwake=rows(conwake);
		nn=rows(nodoswake);
		#Evaluacion de velocidades por metodo explicito de un paso
		#Estado actual
		R=nodoswake(:,1:3);
    #Evaluo FI
    [FI]=FIeval(R,nn,rowwake,conwake,G,niB,njB,nkB,nlB,numelm,Vlibre,wy,deltaT);
    #Convecto todos los nodos de estela
    R=R+FI*deltaT;
    nodoswake=R;
  
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
  	nnodoslibres=length(nodosstartvrtx); #=n+m+1
    #Agrego nuevos nodos del borde libre a la estela
		#Determino posiciones adheridas
		Rsvx=nodos(nodosstartvrtx,1:3);
		#Determino ultimo lugar en la tabla de estela
		nn=rows(nodoswake);
		#Agrego los nodos a la tabla de nodos estela 
		nodoswake(nn+1:nn+nnodoslibres,1:3)=Rsvx;
		#Aumento puntero del tama�o de tabla
		nn=rows(nodoswake);

		#Agrego nuevos anillos de start vortex a tabla de conectividades
		#En la tabla vinculo los nodos de los anillos libres con los nodos del borde libre
    if desprendimiento==1 #Para el caso de convectar solo borde de fuga
      for i=1:n #Fuga
        boundelem=(m-1)*n+i; #anillos adheridos en trailing edge
        wakeelem=rowwake-(n)+i; #anillos de estela en trailing edge
        nodoi=nn-nnodoslibres+i;	#busco directamente los ultimos nodos de nodoswake
        nodoj=nn-nnodoslibres+i+1;
        nodok=conwake(wakeelem,2); #k coincide con el nodo j de wake unido
        nodol=conwake(wakeelem,1); #l coincide con el nodo i de wake unido
        conwake(rowwake+i,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
		  endfor
    elseif desprendimiento==2 #Para el caso de convectar borde de fuga y extremo derecho
      for i=1:m #Derechos
        boundelem=(i-1)*n+n;
        wakeelem=rowwake-(m+n)+i; #vuelvo el puntero desde el final
        nodoi=nn-nnodoslibres+i; #vuelvo el puntero desde el final
        nodol=nn-nnodoslibres+i+1; #vuelvo el puntero desde el final			
        nodoj=conwake(wakeelem,1);#conecta con wake
        nodok=conwake(wakeelem,4);#conecta con wake
        conwake(i+rowwake,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
  		endfor
	  	for i=1:n-1 #Fuga
        boundelem=(m-1)*n+i;
        wakeelem=rowwake-(m+n)+m+i;
        nodoi=nn-nnodoslibres+m+1+i;	
        nodoj=nn-nnodoslibres+m+1+i+1;
        nodok=conwake(wakeelem,2);#conecta con wake
        nodol=conwake(wakeelem,1);#conecta con wake
        conwake(rowwake+m+i,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
  		endfor
	  	i=n; #Esquina
		  boundelem=(m-1)*n+i;
		  wakeelem=rowwake;
		  nodoi=nn-nnodoslibres+(m+1)+n;	
		  nodoj=nn-nnodoslibres+(m+1); #La esquina pertenece al borde derecho
		  nodok=conwake(wakeelem,2);#conecta con wake
		  nodol=conwake(wakeelem,1);#conecta con wake
		  conwake(rowwake+m+i,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
  
    elseif desprendimiento==3 #Borde de fuga y ambos extremos
      r=rowwake+1; #Actualizacion de tabla de conectividad
      for i=1:m #anillos izq
        boundelem=(i-1)*n+1;
        wakeelem=rowwake-(2*m+n)+i; #anillos wake a conectar
        nodoi=conwake(wakeelem,2); #conecta con wake j
        nodoj=nn-nnodoslibres+i; #conecta con adherido
        nodok=nn-nnodoslibres+i+1; #conecta con adherido
        nodol=conwake(wakeelem,3); #conecta con wake k
        conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
        r++;
		  endfor
		for i=1:m #anillos der
      boundelem=i*n;
      wakeelem=rowwake-(m+n)+i; #anillos wake a conectar
			nodoi=nn-nnodoslibres+m+1+i; #conecta con adherido
			nodoj=conwake(wakeelem,1);#conecta con wake i
			nodok=conwake(wakeelem,4);#conecta con wake l
			nodol=nn-nnodoslibres+m+1+i+1; #conecta con adherido
			conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
      r++;
		endfor
    
		for i=1:n-2 #anillos fuga
      boundelem=(m-1)*n+i+1;
      wakeelem=rowwake-(n)+i; #anillos wake a conectar
			nodoi=nn-nnodoslibres+2*(m+1)+i; #conecta con adherido
			nodoj=nn-nnodoslibres+2*(m+1)+i+1; #conecta con adherido
			nodok=conwake(wakeelem,2);#conecta con wake j
			nodol=conwake(wakeelem,1);#conecta con wake i
			conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
      r++;
		endfor

    i=1; #conecto con anillo estela izquierdo
		boundelem=(m-1)*n+i;		
    wakeelemi=rowwake+m; #anillo a la izquierda
    wakeelemd=rowwake+2*m+1; #anillo a la derecha
		nodoi=conwake(wakeelemi,3); #conecto con nodo k
		nodoj=conwake(wakeelemd,1);
		nodok=conwake(wakeelemd,4);
		nodol=conwake(wakeelemi,4); #conecto con nodo l
		conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
    r++;

		i=n; #conecto con anillo estela derecho
		boundelem=(m-1)*n+i;		
    wakeelemd=rowwake+2*m;
    wakeelemi=rowwake+2*m+n-2;
		nodoi=conwake(wakeelemi,2);
		nodoj=conwake(wakeelemd,4); #conecto con nodo l
		nodok=conwake(wakeelemd,3); #conecto con nodo k
		nodol=conwake(wakeelemi,3);
		conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
         
 endif
    


		###TRUNCAMIENTO###############################################################3
    maxconv=130; #numero de tsteps:640 para LC/8 320 para LC/4 160 para LC/2 80 para LC*1 40 para LC*2, 20 para LC*4, 10 para LC*8
		
    if desprendimiento==1		
      nodosxconveccion=(n+1); #nodos convectados en cada conveccion
      anillosxconveccion=n; #anillos convectados en cada conveccion
      if nn >= maxconv*nodosxconveccion #nodos convectados superan maximo numero de nodos?		
        nodoswake(1:end-nodosxconveccion,:)=nodoswake(nodosxconveccion+1:end,:); #Piso el lugar de los nodos lejanos
        nodoswake(end-nodosxconveccion+1:end,:)=[]; #Elimino filas sobrantes
        conwake(1:end- anillosxconveccion,:)=conwake( anillosxconveccion+1:end,:); #Piso el lugar de los anillos lejanos
        conwake(end- anillosxconveccion+1:end,:)=[]; #Elimino filas sobrantes
        conwake(:,1:4)=bsxfun('minus',conwake(:,1:4),nodosxconveccion); #Corrijo numeracion de nodos
      endif

		elseif desprendimiento==2
      nodosxconveccion=((m+1)+n-1); #nodos convectados en cada conveccion
      anillosxconveccion=m+n; #anillos convectados en cada conveccion
      if nn >= maxconv*nodosxconveccion #nodos convectados superan maximo numero de nodos?		
        nodoswake(1:end-nodosxconveccion,:)=nodoswake( nodosxconveccion+1:end,:); #Piso el lugar de los nodos lejanos
        nodoswake(end-nodosxconveccion+1:end,:)=[]; #Elimino filas sobrantes
        conwake(1:end-anillosxconveccion,:)=conwake(anillosxconveccion+1:end,:); #Piso el lugar de los anillos lejanos
        conwake(end-anillosxconveccion+1:end,:)=[]; #Elimino filas sobrantes
        conwake(:,1:4)=bsxfun('minus',conwake(:,1:4),nodosxconveccion); #Corrijo numeracion de nodos
   		endif

    elseif desprendimiento==3
      nodosxconveccion=(2*(m+1)+n-2); #nodos convectados en cada conveccion
      anillosxconveccion=2*m+n; #anillos convectados en cada conveccion
       if nn >= maxconv*nodosxconveccion #nodos convectados superan maximo numero de nodos?		
        nodoswake(1:end-nodosxconveccion,:)=nodoswake(nodosxconveccion+1:end,:); #Piso el lugar de los nodos lejanos
        nodoswake(end-nodosxconveccion+1:end,:)=[]; #Elimino filas sobrantes
        conwake(1:end-anillosxconveccion,:)=conwake(anillosxconveccion+1:end,:); #Piso el lugar de los anillos lejanos
        conwake(end-anillosxconveccion+1:end,:)=[]; #Elimino filas sobrantes
        conwake(:,1:4)=bsxfun('minus',conwake(:,1:4),nodosxconveccion); #Corrijo numeracion de nodos
   		endif
		endif
    
	endif

################################################################################################		
################################################################################################
	if initestela==1 #No hay estela, paso inicial, convecto nodos libres a estela
		
  	q=1; #Inicio puntero que cuenta nodos generados en estela
		
    if desprendimiento==1 #Para el caso de convectar solo borde de fuga
      #Busco nodos del borde libre
      nodosbordefuga=[m*(n+1)+1:(m+1)*(n+1)];#con esquinas
		  nodosstartvrtx=[nodosbordefuga];
		  nnodoslibres=size(nodosstartvrtx,2);
      #Convecto los nodos en el borde libre
	    for i=nodosstartvrtx #Calculo de velocidades
        #Asigno variables
        R=nodos(i,1:3);
        V=AeroVinducida_rev5(R,G,niB,njB,nkB,nlB,1);	#Velocidad inducida por superficie
        #Agrego velocidad libre
        V=V+Vlibre;
        #Expreso en sistema B
        k1=-wy*[R(1,3),0,-R(1,1)];
        Raux=R+0.5*k1*deltaT;
        k2=-wy*[Raux(1,3),0,-Raux(1,1)];
        Raux=R+0.5*k2*deltaT;
        k3=-wy*[Raux(1,3),0,-Raux(1,1)];
        Raux=R+k3*deltaT;
        k4=-wy*[Raux(1,3),0,-Raux(1,1)];
        V=V+(k1+2*k2+2*k3+k4)/6;
        #Asigno en tabla, nodos convectados segun orden nodosstartvrtx
        nodoswake(q,1:3)=[R+V*deltaT]; #Primero nodos convectados
        nodoswake(q+nnodoslibres,1:3)=[R]; #Despues agrego nodos todavia adheridos
        q++;
      endfor    
      for i=1:n #Actualizacion de tabla de conectividad
        boundelem=(m-1)*n+i; #ultima fila
        wakeelem=i; #no hay elementos previos en estela
        #Los nodos en tabla nodoswake est�n ordenados seg�n nodosstartvrtx. Primero se convectan los kl y luego se agregan los ij.
        nodosconvectados=n+1; #solo borde de fuga
        nodoi=nodosconvectados+i;
        nodoj=nodosconvectados+i+1;
        nodok=i+1; 
        nodol=i;
        conwake(i,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
      endfor
    
    elseif desprendimiento==2      
      #Busco nodos del borde libre
		  nodosbordeder=[n+1:n+1:(m+1)*(n+1)];
		  nodosbordefuga=[m*(n+1)+1:(m+1)*(n+1)-1];#sin esquinas
		  nodosstartvrtx=[nodosbordeder,nodosbordefuga];
		  nnodoslibres=size(nodosstartvrtx,2);
      #Convecto los nodos en el borde libre
	    for i=nodosstartvrtx #Calculo de velocidades
        R=nodos(i,1:3);
        V=AeroVinducida_rev5(R,G,niB,njB,nkB,nlB,1);	#Velocidad inducida por superficie
        #Agrego velocidad libre
        V=V+Vlibre;
        #Expreso en sistema B
        k1=-wy*[R(1,3),0,-R(1,1)];
        Raux=R+0.5*k1*deltaT;
        k2=-wy*[Raux(1,3),0,-Raux(1,1)];
        Raux=R+0.5*k2*deltaT;
        k3=-wy*[Raux(1,3),0,-Raux(1,1)];
        Raux=R+k3*deltaT;
        k4=-wy*[Raux(1,3),0,-Raux(1,1)];
        V=V+(k1+2*k2+2*k3+k4)/6;
        #Asigno en tabla
        nodoswake(q,1:3)=[R+V*deltaT]; #Primero nodos convectados
        nodoswake(q+nnodoslibres,1:3)=[R]; #Despues agrego nodos todavia adheridos
        q++;
  		endfor    
      #2 - armo matriz conectividad 
      r=1; #Actualizacion de tabla de conectividad
      for q=1:m #elementos der
        nodoi=nnodoslibres+q;
        nodoj=q;
        nodok=q+1;
        nodol=nnodoslibres+q+1;
        conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(q*n)];
        r++;
      endfor
      for q=1:n-1 #elementos fuga
        boundelem=(m-1)*n+q;
        nodoi=nnodoslibres+m+1+q;
        nodoj=nnodoslibres+m+1+q+1;
        nodok=m+1+q+1;
        nodol=m+1+q;
        conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
        r++;  
      endfor
      q=n; #equivalencio con anillo estela derecho
      boundelem=(m-1)*n+q;			
      nodoi=nnodoslibres+m+1+q;
      nodoj=conwake(m,4); #equivalencio
      nodok=conwake(m,3); #equivalencio
      nodol=m+1+q;
      conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
      
  elseif desprendimiento==3
      #Busco los nodos adheridos en el borde libre
      nodosbordeizq=[1:n+1:m*(n+1)+1];
      nodosbordeder=[n+1:n+1:(m+1)*(n+1)];
		  nodosbordefuga=[m*(n+1)+2:(m+1)*(n+1)-1];#sin esquinas
		  nodosstartvrtx=[nodosbordeizq,nodosbordeder,nodosbordefuga];
		  nnodoslibres=size(nodosstartvrtx,2);
      #Convecto los nodos en el borde libre
	    for i=nodosstartvrtx #Calculo de velocidades
        R=nodos(i,1:3);
        V=AeroVinducida_rev5(R,G,niB,njB,nkB,nlB,1);	#Velocidad inducida por superficie
        #Agrego velocidad libre
        V=V+Vlibre;
        #Expreso en sistema B
        k1=-wy*[R(1,3),0,-R(1,1)];
        Raux=R+0.5*k1*deltaT;
        k2=-wy*[Raux(1,3),0,-Raux(1,1)];
        Raux=R+0.5*k2*deltaT;
        k3=-wy*[Raux(1,3),0,-Raux(1,1)];
        Raux=R+k3*deltaT;
        k4=-wy*[Raux(1,3),0,-Raux(1,1)];
        V=V+(k1+2*k2+2*k3+k4)/6;
        #Asigno en tabla
        nodoswake(q,1:3)=[R+V*deltaT]; #Primero nodos convectados
        nodoswake(q+nnodoslibres,1:3)=[R]; #Despues agrego nodos todavia adheridos
        q++;
  		endfor
      #2 - armo matriz conectividad 
      r=1; #Actualizacion de tabla de conectividad
      for q=1:m #anillos izq
        boundelem=(q-1)*n+1;
        nodoi=q;
        nodoj=nnodoslibres+q;
        nodok=nnodoslibres+q+1;
        nodol=q+1;
        conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
        r++;
      endfor
      for q=1:m #anillos der
        boundelem=q*n;
        nodoi=(nnodoslibres+m+1)+q;
        nodoj=(m+1)+q;
        nodok=(m+1)+q+1;
        nodol=(nnodoslibres+m+1)+q+1;
        conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
        r++;
      endfor
      for q=2:n-1 #anillos fuga
        boundelem=(m-1)*n+q-1;
        nodoi=nnodoslibres+2*(m+1)+q-1;
        nodoj=nnodoslibres+2*(m+1)+q;
        nodok=2*(m+1)+q;
        nodol=2*(m+1)+q-1;
        conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
        r++;  
      endfor
      q=1; #conecto con anillo estela izquierdo
      boundelem=(m-1)*n+q		
      nodoi=conwake(m,3) #conecto con nodo k
      nodoj=conwake(2*m+1,1);
      nodok=conwake(2*m+1,4);
      nodol=conwake(m,4) #conecto con nodo l
      conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
      r++;
      q=n #conecto con anillo estela derecho
      boundelem=(m-1)*n+q;		
      nodoi=conwake(2*m+n-2,2);
      nodoj=conwake(2*m,4); #conecto con nodo l
      nodok=conwake(2*m,3); #conecto con nodo k
      nodol=conwake(2*m+n-2,3);
      conwake(r,:)=[nodoi,nodoj,nodok,nodol,G(boundelem)];
    endif
endif
    
endfunction
