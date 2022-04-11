function [u,up,upp,v,vp,vpp,w,wp,wpp,t,tp,tpp,nax,nflap,nedge,ntor]=modoseval_r1(x,L,nax,nflap)
## Evaluacion de funciones de forma modales (Libro de Hodges)

#Modos de estiramiento de barra fija-libre
for i=1:nax
	u(i,:)=[sin( (i-0.5)*pi*x/L)];
	up(i,:)=[(i-0.5)*pi/L*cos( (i-0.5)*pi*x/L) ];
	upp(i,:)=[-((i-0.5)*pi/L)^2*sin( (i-0.5)*pi*x/L)];
endfor

ntor=nax;
t=u;
tp=up;
tpp=upp;

#Modos de flexion de viga empotrada-libre 
alfaiL=1.875104068711961166;#1.8751
alfai=alfaiL/L;
betai=(cosh(alfaiL)+cos(alfaiL)) / ( sinh(alfaiL)+sin(alfaiL) ); #betai=0.734096;
v(1,:)=[cosh(alfai*x)-cos(alfai*x)-betai*(sinh(alfai*x)-sin(alfai*x) )];
vp(1,:)=alfai*[sinh(alfai*x)+sin(alfai*x)-betai*(cosh(alfai*x)-cos(alfai*x) )];
vpp(1,:)=alfai^2*[cosh(alfai*x)+cos(alfai*x)-betai*(sinh(alfai*x)+sin(alfai*x) )];

if nflap>1
	alfaiL=4.694091132974174576;#4.69409
	alfai=alfaiL/L;	#betai=1.01847;
	betai=(cosh(alfaiL)+cos(alfaiL)) / ( sinh(alfaiL)+sin(alfaiL) );
	v(2,:)=[cosh(alfai*x)-cos(alfai*x)-betai*(sinh(alfai*x)-sin(alfai*x) )];
	vp(2,:)=alfai*[sinh(alfai*x)+sin(alfai*x)-betai*(cosh(alfai*x)-cos(alfai*x) )];
	vpp(2,:)=alfai^2*[cosh(alfai*x)+cos(alfai*x)-betai*(sinh(alfai*x)+sin(alfai*x) )];
endif
if nflap>2
	alfaiL=7.854757438237612565;#7.85476
	alfai=alfaiL/L;	#betai=0.999224
	betai=(cosh(alfaiL)+cos(alfaiL)) / ( sinh(alfaiL)+sin(alfaiL) );
	v(3,:)=[cosh(alfai*x)-cos(alfai*x)-betai*(sinh(alfai*x)-sin(alfai*x) )];
	vp(3,:)=alfai*[sinh(alfai*x)+sin(alfai*x)-betai*(cosh(alfai*x)-cos(alfai*x) )];
	vpp(3,:)=alfai^2*[cosh(alfai*x)+cos(alfai*x)-betai*(sinh(alfai*x)+sin(alfai*x) )];
endif
if nflap>3
	alfaiL=10.99554073487546699;#7.85476
	alfai=alfaiL/L;	#betai=0.999224
	betai=(cosh(alfaiL)+cos(alfaiL)) / ( sinh(alfaiL)+sin(alfaiL) );
	v(4,:)=[cosh(alfai*x)-cos(alfai*x)-betai*(sinh(alfai*x)-sin(alfai*x) )];
	vp(4,:)=alfai*[sinh(alfai*x)+sin(alfai*x)-betai*(cosh(alfai*x)-cos(alfai*x) )];
	vpp(4,:)=alfai^2*[cosh(alfai*x)+cos(alfai*x)-betai*(sinh(alfai*x)+sin(alfai*x) )];
endif
if nflap>4
	alfaiL=14.13716839104647058;#7.85476
	alfai=alfaiL/L;	#betai=0.999224
	betai=(cosh(alfaiL)+cos(alfaiL)) / ( sinh(alfaiL)+sin(alfaiL) );
	v(5,:)=[cosh(alfai*x)-cos(alfai*x)-betai*(sinh(alfai*x)-sin(alfai*x) )];
	vp(5,:)=alfai*[sinh(alfai*x)+sin(alfai*x)-betai*(cosh(alfai*x)-cos(alfai*x) )];
	vpp(5,:)=alfai^2*[cosh(alfai*x)+cos(alfai*x)-betai*(sinh(alfai*x)+sin(alfai*x) )];
endif

if nflap>5
	alfaiL=17.27875953208823633;#7.85476
	alfai=alfaiL/L;	#betai=0.999224
	betai=(cosh(alfaiL)+cos(alfaiL)) / ( sinh(alfaiL)+sin(alfaiL) );
	v(6,:)=[cosh(alfai*x)-cos(alfai*x)-betai*(sinh(alfai*x)-sin(alfai*x) )];
	vp(6,:)=alfai*[sinh(alfai*x)+sin(alfai*x)-betai*(cosh(alfai*x)-cos(alfai*x) )];
	vpp(6,:)=alfai^2*[cosh(alfai*x)+cos(alfai*x)-betai*(sinh(alfai*x)+sin(alfai*x) )];
endif

if nflap>6
	for i=7:nflap
		alfaiL=(2*i-1)*pi/2;
		alfai=alfaiL/L;
		betai=(cosh(alfaiL)+cos(alfaiL)) / ( sinh(alfaiL)+sin(alfaiL) );
		v(i,:)=[cosh(alfai*x)-cos(alfai*x)-betai*(sinh(alfai*x)-sin(alfai*x) )];
		vp(i,:)=alfai*[sinh(alfai*x)+sin(alfai*x)-betai*(cosh(alfai*x)-cos(alfai*x) )];
		vpp(i,:)=alfai^2*[cosh(alfai*x)+cos(alfai*x)-betai*(sinh(alfai*x)+sin(alfai*x) )];
	endfor
endif

#Reemplazo para puente con modos apoyado-apoyado
#for i=1:nflap
#	v(i,:)=sin(i*pi/L*x);
#	vp(i,:)=(i*pi/L)*cos(i*pi/L*x);
#	vpp(i,:)=-(i*pi/L)^2*sin(i*pi/L*x);
#endfor

nedge=nflap;
w=v;
wp=vp;
wpp=vpp;



endfunction

