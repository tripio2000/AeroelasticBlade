function [V]=AeroDblPnt(G,Rcp,nhat,A,Reval)
## Calcula velocidades inducidas por doblete puntual
#Rcp punto de control del panel
#nhat normal al panel
#A es el area del panel
#Reval punto a evaluar velocidad
#G circulacion del anillo
#V es la velocidad inducida por doblete

%	#Preasigno tama�os		
%		V=zeros(size(G,1),3); 
%		r1=V;
%		nr1=V;
%		r2=V;
%		nr2=V;
%		L=V;
%		k=V;
%		nk=V;
%		e12=V;
    		
x0=Rcp(:,1);
y0=Rcp(:,2);
z0=Rcp(:,3);
nx=nhat(:,1);
ny=nhat(:,2);
nz=nhat(:,3);
x=Reval(:,1);
y=Reval(:,2);
z=Reval(:,3);
r=Reval-Rcp;
nr=norm(r,'rows');

d=(nx.*r(:,1)+ny.*r(:,2)+nz.*r(:,3));

%u=nx.*nr.^2-3*r(:,1).*d;
%v=ny.*nr.^2-3*r(:,2).*d;
%w=nz.*nr.^2-3*r(:,3).*d;
%V=[u,v,w].*(-G./A)./(4*pi*nr.^5);

nr(nr<1E-15)=inf;

u=nx./(nr.^3)-3*r(:,1).*d./(nr.^5);
v=ny./(nr.^3)-3*r(:,2).*d./(nr.^5);
w=nz./(nr.^3)-3*r(:,3).*d./(nr.^5);
V=[u,v,w].*(-G./(4*pi*A));


endfunction
