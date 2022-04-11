# Rotación propia de vectores
# q es el vector con los ángulos de rotación (x,y,z)
#Los vectores de entrada y salida son fila (horizontales)

function [vectorout]=xformNB(q,vectorin)
	mx=cos(q(1));
	nx=sin(q(1));
	R1=[1,0,0;0,mx,-nx;0,nx,mx];
	my=cos(q(2));
	ny=sin(q(2));
	R2=[my,0,ny;0,1,0;-ny,0,my];
	mz=cos(q(3));
	nz=sin(q(3));
	R3=[mz,-nz,0;nz,mz,0;0,0,1];
	vectorout=(R3*R2*R1*vectorin');
		# vectorout=(R3*R2*R1*[vectorin(:,1);vectorin(:,2);vectorin(:,3)]);
		vectorout=vectorout';
endfunction
