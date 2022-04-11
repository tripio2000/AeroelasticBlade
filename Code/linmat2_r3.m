function [M,Kc,Fc,Ipolar,Ke]=linmat2_r3(h,x1,x2,deltax,BMassinterp1,Izzinterp1,Iyyinterp1,Iyzinterp1,Syinterp1,Szinterp1,EAinterp1,EJzzinterp1,EJyzinterp1,EJyyinterp1,GJpinterp1,BMassinterp2,Izzinterp2,Iyyinterp2,Iyzinterp2,Syinterp2,Szinterp2,EAinterp2,EJzzinterp2,EJyzinterp2,EJyyinterp2,GJpinterp2,ui1,ui2,vi1,vi2,wi1,wi2,ti1,ti2,uip1,uip2,vip1,vip2,wip1,wip2,vipp1,vipp2,wipp1,wipp2,tip1,tip2,nax,nflap,ntor);
## Evaluo las matrices lineales M,G,Kc,Ke,Fc,Fg,Ipolar
## uso: [M,Kc,Fc,Ipolar,Ke]=linmat2_r1(h,x1,x2,deltax,BMassinterp1,Izzinterp1,Iyyinterp1,Iyzinterp1,Syinterp1,Szinterp1,EAinterp1,EJzzinterp1,EJyzinterp1,EJyyinterp1,GJpinterp1,BMassinterp2,Izzinterp2,Iyyinterp2,Iyzinterp2,Syinterp2,Szinterp2,EAinterp2,EJzzinterp2,EJyzinterp2,EJyyinterp2,GJpinterp2,ui1,ui2,vi1,vi2,wi1,wi2,ti1,ti2,uip1,uip2,vip1,vip2,wip1,wip2,vipp1,vipp2,wipp1,wipp2,tip1,tip2);
## Agrego momentos estaticos inerciales
%% Corrijo error de parentesis en masa torsion
## Rev3: corrección en teoría, cambio matrices M,Kc
## Author: NGT

%nax=size(ui1,1);
%nflap=size(vi1,1);
%ntor=size(ti1,1);

#Evaluo la matriz M
      M=zeros(nax+2*nflap+ntor,nax+2*nflap+ntor);
      #Muu
      for i=1:nax
        for j=1:nax
          M(i,j)=deltax/2*sum(BMassinterp1.*ui1(i,:).*ui1(j,:)+BMassinterp2.*ui2(i,:).*ui2(j,:)); #m*fiuu
        endfor
      endfor
      #Muv,Muw
      for i=1:nax
        for j=1:nflap
          M(i,j+nax)=-deltax/2*sum(Szinterp1.*ui1(i,:).*vip1(j,:)+Szinterp2.*ui2(i,:).*vip2(j,:)); #-Sz*fiuv'
          M(i,j+nflap+nax)=-deltax/2*sum(Syinterp1.*ui1(i,:).*wip1(j,:)+Syinterp2.*ui2(i,:).*wip2(j,:)); #-Sy*fiuw'
        endfor
      endfor  
      M(nax+1:nax+2*nflap,1:nax)=M(1:nax,nax+1:nax+2*nflap)'; #[Mvu;Mwu]           
      #Mvv,Mww
      for i=1:nflap
        for j=1:nflap
          M(i+nax,j+nax)=deltax/2*sum(BMassinterp1.*vi1(i,:).*vi1(j,:)+Izzinterp1.*vip1(i,:).*vip1(j,:)+BMassinterp2.*vi2(i,:).*vi2(j,:)+Izzinterp2.*vip2(i,:).*vip2(j,:)); #m*fivv+Izz*fiv'v'
          M(i+nflap+nax,j+nflap+nax)=deltax/2*sum(BMassinterp1.*wi1(i,:).*wi1(j,:)+Iyyinterp1.*wip1(i,:).*wip1(j,:)+BMassinterp2.*wi2(i,:).*wi2(j,:)+Iyyinterp2.*wip2(i,:).*wip2(j,:)); #m*fiww+Iyy*fiw'w'
        endfor
      endfor
      #Mvw
      for i=1:nflap
        for j=1:nflap
          M(i+nax,j+nflap+nax)=deltax/2*sum(Iyzinterp1.*vip1(i,:).*wip1(j,:)+Iyzinterp2.*vip2(i,:).*wip2(j,:)); #Iyz*fiv'w'
        endfor
      endfor
      M(nflap+nax+1:nax+2*nflap,nax+1:nax+nflap)=M(nax+1:nax+nflap,nflap+nax+1:nax+2*nflap);
      #Mvt,Mwt
      for i=1:nflap
        for j=1:ntor
          M(nax+i,nax+2*nflap+j)=-deltax/2*sum(Syinterp1.*vi1(i,:).*ti1(j,:)+Syinterp2.*vi2(i,:).*ti2(j,:)); #-Sy*fivt
          M(nax+nflap+i,nax+2*nflap+j)=deltax/2*sum(Szinterp1.*wi1(i,:).*ti1(j,:)+Szinterp2.*wi2(i,:).*ti2(j,:)); #+Sz*fiwt
        endfor
      endfor
      M(nax+2*nflap+1:end,nax+1:nax+2*nflap)=M(nax+1:nax+2*nflap,nax+2*nflap+1:end)'; #[Mtv;Mtw]
	    #Mtt
      for i=1:ntor
        for j=1:ntor
          M(nax+2*nflap+i,nax+2*nflap+j)=deltax/2*sum((Iyyinterp1+Izzinterp1).*ti1(i,:).*ti1(j,:)+(Iyyinterp2+Izzinterp2).*ti2(i,:).*ti2(j,:)); #(Iyy+Izz)*fitt
        endfor
      endfor
      
      
	#Evaluo la matriz Kc
      Kc=zeros(nax+2*nflap+ntor,nax+2*nflap+ntor);
%      #Axial-axial
%      for i=1:nax
%        for j=1:nax
%          Kc(i,j)=deltax/2*sum(BMassinterp1.*ui1(i,:).*ui1(j,:)+BMassinterp2.*ui2(i,:).*ui2(j,:)); #Kcuu
%        endfor
%      endfor
%      #Torsion-torsion
%      for i=1:ntor
%        for j=1:ntor
%          Kc(nax+2*nflap+i,nax+2*nflap+j)=deltax/2*sum(Izzinterp1.*ti1(i,:).*ti1(j,:)+Izzinterp2.*ti2(i,:).*ti2(j,:));#Kctt
%        endfor
%      endfor
%      
%      #Flexion-flexion
%      for i=1:nflap
%        for j=1:nflap
%          Kc(nax+i,nax+j)=deltax/2*sum(Izzinterp1.*vip1(i,:).*vip1(j,:)+Izzinterp2.*vip2(i,:).*vip2(j,:)); #Kcvv
%          Kc(nax+i,nax+nflap+j)=deltax/2*sum(Iyzinterp1.*vip1(i,:).*wip1(j,:)+Iyzinterp2.*vip2(i,:).*wip2(j,:));	#Kcvw	
%          Kc(nax+nflap+i,nax+nflap+j)=deltax/2*sum(BMassinterp1.*wi1(i,:).*wi1(j,:)+BMassinterp2.*wi2(i,:).*wi2(j,:)+Iyyinterp1.*wip1(i,:).*wip1(j,:)+Iyyinterp2.*wip2(i,:).*wip2(j,:)); #Kcww
%        endfor
%      endfor
%      Kc(nax+nflap:end-1,nax+1:nax+nflap)=Kc(nax+1:nax+nflap,nax+nflap:end-1)';#Kcwv
%
%      #Axial-flap
%      for i=1:nax
%        for j=1:nflap
%          Kc(i,nax+j)=deltax/2*sum(Szinterp1.*ui1(i,:).*vip1(j,:)+Szinterp2.*ui2(i,:).*vip2(j,:)); #Kcuv
%        endfor
%      endfor
%      Kc(nax+1:nax+nflap,1:nax)=Kc(1:nax,nax+1:nax+nflap)'; #Kcvu
%      
%      #Flexion-torsion
%      for i=1:nflap
%        for j=1:ntor
%          Kc(nax+i,nax+2*nflap+j)=deltax/2*sum(Syinterp1.*vi1(i,:).*ti1(j,:)+Syinterp2.*vi2(i,:).*ti2(j,:)); #Kcvt
%          Kc(nax+nflap+i,nax+2*nflap+j)=-deltax/2*sum(Szinterp1.*wi1(i,:).*ti1(j,:)+Szinterp2.*wi2(i,:).*ti2(j,:)); #Kcwt
%        endfor
%      endfor
%      Kc(nax+2*nflap+1:end,nax+1:nax+2*nflap)=Kc(nax+1:nax+2*nflap,nax+2*nflap+1:end)'; #[Kctv;Kctw]
 
     #Kcuu
      for i=1:nax
        for j=1:nax
          Kc(i,j)=deltax/2*sum(BMassinterp1.*ui1(i,:).*ui1(j,:)+BMassinterp2.*ui2(i,:).*ui2(j,:)); #m*fiuu
        endfor
      endfor
      #Kcuv,Kcuw
      for i=1:nax
        for j=1:nflap
          Kc(i,j+nax)=-deltax/2*sum(Szinterp1.*ui1(i,:).*vip1(j,:)+Szinterp2.*ui2(i,:).*vip2(j,:)); #-Sz*fiuv'
          Kc(i,j+nflap+nax)=-deltax/2*sum(Syinterp1.*ui1(i,:).*wip1(j,:)+Syinterp2.*ui2(i,:).*wip2(j,:)); #-Sy*fiuw'
        endfor
      endfor  
      Kc(nax+1:nax+2*nflap,1:nax)=Kc(1:nax,nax+1:nax+2*nflap)'; #[Mvu;Mwu]           
      #Kcww
      for i=1:nflap
        for j=1:nflap
          Kc(i+nflap+nax,j+nflap+nax)=deltax/2*sum(BMassinterp1.*wi1(i,:).*wi1(j,:)+BMassinterp2.*wi2(i,:).*wi2(j,:)); #m*fiww+Iyy*fiw'w'
        endfor
      endfor
      #Kcwt
      for i=1:nflap
        for j=1:ntor
          Kc(nax+nflap+i,nax+2*nflap+j)=deltax/2*sum(Szinterp1.*wi1(i,:).*ti1(j,:)+Szinterp2.*wi2(i,:).*ti2(j,:)); #+Sz*fiwt
        endfor
      endfor
      Kc(nax+2*nflap+1:end,nax+1:nax+2*nflap)=Kc(nax+1:nax+2*nflap,nax+2*nflap+1:end)'; #[Mtv;Mtw]
	    
  #Evaluo el vector fc
      Fc=zeros(nax+2*nflap+ntor,1);
      for i=1:nax
        Fc(i,1)=deltax/2*sum(BMassinterp1.*(x1+h(ones(size(x1)))).*ui1(i,:)+BMassinterp2.*(x2+h(ones(size(x2)))).*ui2(i,:));
      endfor	
      for i=1:nflap
        Fc(nax+i,1)=-deltax/2*sum(Szinterp1.*(x1+h(ones(size(x1)))).*vip1(i,:)+Szinterp2.*(x2+h(ones(size(x2)))).*vip2(i,:));
    %		Fc(nax+nflap+i,1)=-deltax/2*sum(Syinterp1.*(x1+h(ones(size(x1)))).*wip1(i,:)-wi1(i,:)+Syinterp2.*(x2+h(ones(size(x2)))).*wip2(i,:)-wi2(i,:)); %Revisado parentesis!!!!
        Fc(nax+nflap+i,1)=-deltax/2*sum(Syinterp1.*((x1+h(ones(size(x1)))).*wip1(i,:)-wi1(i,:))+Syinterp2.*((x2+h(ones(size(x2)))).*wip2(i,:)-wi2(i,:)));
      endfor	
      for i=1:ntor
        Fc(nax+2*nflap+i,1)=deltax/2*sum(Iyzinterp1.*ti1(i,:)+Iyzinterp2.*ti2(i,:));
      endfor	
	
	Ipolar=deltax/2*sum((0.5*Iyyinterp1+0.5*Iyyinterp2)*2+(0.5*BMassinterp1+0.5*BMassinterp2).*(x1.^2+x2.^2));

#Evaluo la matriz el�stica lineal
      Ke=zeros(nax+2*nflap+ntor,nax+2*nflap+ntor);
      for i=1:nax
        for j=1:nax
          Ke(i,j)=sum(deltax/2*(EAinterp1.*uip1(i,:).*uip1(j,:)+EAinterp2.*uip2(i,:).*uip2(j,:))); #Keuu
        endfor
      endfor

      for i=1:ntor
        for j=1:ntor
          Ke(nax+2*nflap+i,nax+2*nflap+j)=sum(deltax/2*(GJpinterp1.*tip1(i,:).*tip1(j,:)+GJpinterp2.*tip2(i,:).*tip2(j,:))); #Kett
        endfor
      endfor
      
      for i=1:nflap
        for j=1:nflap
          Ke(i+nax,j+nax)=sum(deltax/2*(EJzzinterp1.*vipp1(i,:).*vipp1(j,:)+EJzzinterp2.*vipp2(i,:).*vipp2(j,:))); #Kevv
          Ke(i+nax,j+nflap+nax)=sum(deltax/2*(EJyzinterp1.*vipp1(i,:).*wipp1(j,:)+EJyzinterp2.*vipp2(i,:).*wipp2(j,:))); #Kevw		
          Ke(i+nflap+nax,j+nflap+nax)=sum(deltax/2*(EJyyinterp1.*wipp1(i,:).*wipp1(j,:)+EJyyinterp2.*wipp2(i,:).*wipp2(j,:))); #Keww		
        endfor
      endfor
      Ke(nax+nflap+1:nax+2*nflap,nax+1:nax+nflap)=Ke(nax+1:nax+nflap,nax+nflap+1:nax+2*nflap)'; #Kewv


endfunction