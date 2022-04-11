function [V]=AeroVorSegGPU(G,Ra,Rb,Rc)
## Evaluation of velocity induced by a single vortex segment.
## Ra,Rb are the 3D position vectors of the segment ends.
## Rc is the 3D position vector of the evaluation point.
## G is the vortex intensity.
## V is the induced velocity.
## Ra,Rb,Rc have dimensions [n,3]
## G has dimension [n,1]
## function call: [V]=vorseg_gpu_rev1 (G,Ra,Rb,Rc)

    
# Transfer mat to OpenCL memory
Ra = oclArray (Ra);
Rb = oclArray (Rb);
Rc = oclArray (Rc);
G = oclArray (G);
   		
    r1=Rc-Ra;
    r2=Rc-Rb;
    L=Rb-Ra;
    
%    nr1=norm(r1,2,'rows');
%    nr2=norm(r2,2,'rows');
%    nL=norm(L,2,'rows');
%    k=cross(L,r1,2);
%		 nk=norm(k,2,'rows');

    nr1=sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
    nr2=sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2);
    nL=( L(:,1).^2+L(:,2).^2+L(:,3).^2);		
    #k=[L(:,2).*r1(:,3)-r1(:,2).*L(:,3), -(L(:,1).*r1(:,3)-r1(:,1).*L(:,3)),L(:,1).*r1(:,2)-r1(:,1).*L(:,2)];
    k1=L(:,2).*r1(:,3)-r1(:,2).*L(:,3);
    k2=-(L(:,1).*r1(:,3)-r1(:,1).*L(:,3));
    k3=L(:,1).*r1(:,2)-r1(:,1).*L(:,2);
    nk=( k1.^2+k2.^2+k3.^2 );		
    
    vindmax=inf;
    
    cutoff=abs(G)/(2*pi*vindmax);
    nk=nk+(cutoff).^2;
    
	  nk(nk<1E-15)=inf;   #when L and r1 are parallel
		nr1(nr1<1E-15)=inf; #when evaluation point is close to segment extreme a
		nr2(nr2<1E-15)=inf; #when evaluation point is close to segment extreme b

##		e12=(r1./nr1-r2./nr2);
    e12_1=(r1(:,1)./nr1-r2(:,1)./nr2);
    e12_2=(r1(:,2)./nr1-r2(:,2)./nr2);
    e12_3=(r1(:,3)./nr1-r2(:,3)./nr2);
    
    dotprodLe12=L(:,1).*e12_1+L(:,2).*e12_2+L(:,3).*e12_3;
    
##		V=1/(4*pi)*G.*dotprodLe12.*(k./nk);
	  V1=1/(4*pi)*G.*dotprodLe12.*(k1./nk);
    V2=1/(4*pi)*G.*dotprodLe12.*(k2./nk);
    V3=1/(4*pi)*G.*dotprodLe12.*(k3./nk);
  
    # transfer results back to octave memory
    V1 = ocl_to_octave (V1);
    V2 = ocl_to_octave (V2);
    V3 = ocl_to_octave (V3);
    V=[V1,V2,V3];
    
    #Clear GPU memmory
    clear -x V

endfunction
