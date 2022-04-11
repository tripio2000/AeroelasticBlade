function [V]=AeroVorSeg_rev1 (G,Ra,Rb,Rc)
## funcion vectorizada
## C�lculo de velocidad inducida por segmento vorticoso
## Ra,Rb son los vectores posicion que definen el segmento. 
## Rc es el vector posicion de evaluacion
## G es la intensidad del vortice
## 16-3-15 armo listas de campo cercano y lejano
## rev1: 5/12/18 recibe tamaño del anillo "a" para usar el mismo cutoff en todo el anillo
# function [V]=vorseg_bsx_rev1 (G,Ra,Rb,Rc,La)
    
		#Preasigno tama�os		
		V=zeros(size(G,1),3); 
		r1=V;
		nr1=V;
		r2=V;
		nr2=V;
		L=V;
		k=V;
		nk=V;
		e12=V;

    		
%    r1=Rc-Ra;
%    nr1=norm(r1,2,'rows');
%    r2=Rc-Rb;
%    nr2=norm(r2,2,'rows');
%    L=Rb-Ra;
%    nL=norm(L,2,'rows');
%    k=cross(L,r1,2);
%		nk=norm(k,2,'rows');

		r1=bsxfun(@minus,Rc,Ra);
		nr1=sqrt( bsxfun(@plus, bsxfun(@plus,r1(:,1).^2,r1(:,2).^2),r1(:,3).^2 ) );
		r2=bsxfun(@minus,Rc,Rb);
		nr2=sqrt( bsxfun(@plus, bsxfun(@plus,r2(:,1).^2,r2(:,2).^2),r2(:,3).^2 ) );		
		L=bsxfun(@minus,Rb,Ra);
		nL=( bsxfun(@plus, bsxfun(@plus,L(:,1).^2,L(:,2).^2),L(:,3).^2 ) );		
		k=[(L(:,2).*r1(:,3)-r1(:,2).*L(:,3)), -(L(:,1).*r1(:,3)-r1(:,1).*L(:,3)),(L(:,1).*r1(:,2)-r1(:,1).*L(:,2))];
		nk=( bsxfun(@plus, bsxfun(@plus,k(:,1).^2,k(:,2).^2),k(:,3).^2 ) );		

%    cutoff=1E-1; %paper leishman y Bagai, de obs fisica el radio del vortice es proximo al espesor del perfil. 
    #potencial, 1E-3,1E-2

    %    vindmax=102/10; #(0.3*a)/LC para que sea incompresible, 1 para que sea = vref
%    vindmax=102/47.744; #wy=2.14;
%    vindmax=102/95.488; #wy=2*2.14;
%    vindmax=102/119.36; #wy=2.5*2.14;
%    vindmax=102/42.455; #wy=2;
     vindmax=inf;
    
    cutoff=abs(G)/(2*pi*vindmax);
    
    nk=nk+(cutoff).^2;
    
%    nk=nk+(cutoff.*nL).^2;
%		nk=nk+(cutoff*La).^2;	
%	  nk(nk<cutoff*nL)=inf; #cuando L y r1 son paralelos

	  nk(nk<1E-15)=inf; #cuando L y r1 son paralelos
		nr1(nr1<1E-15)=inf; #cuando el punto está próximo al extremo a
		nr2(nr2<1E-15)=inf; #cuando el punto está próximo al extremo b

		e12=(r1./nr1-r2./nr2);
    
		V=1/(4*pi)*G.*(k./nk).*dot(L,e12,2);
	
endfunction
