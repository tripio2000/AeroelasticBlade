function [FI]=FIeval(R,nn,rowwake,conwake,G,niB,njB,nkB,nlB,numelm,Vlibre,wy,deltaT)
   knull=zeros(nn,1);
   wiB=R(conwake(:,1),:);
	 wjB=R(conwake(:,2),:);
	 wkB=R(conwake(:,3),:);
	 wlB=R(conwake(:,4),:);
	 Gw=conwake(:,5);

   method=1; #Euler simple 1, RK2 2, RK4 4
    V1=zeros(nn,3);
if method==1   
    V1=V1+AeroVinducida_rev5(R,G,niB,njB,nkB,nlB,1);	#Velocidad inducida por superficie
    V1=V1+AeroVinducida_rev5(R,Gw,wiB,wjB,wkB,wlB,1); #Velocidad inducida por estela
    V1=V1+repmat(Vlibre,nn,1); #Velocidad libre y rotacion
    k1=-wy*[R(:,3),knull,-R(:,1)];
		Raux=R+0.5*k1*deltaT;
		k2=-wy*[Raux(:,3),knull,-Raux(:,1)];
		Raux=R+0.5*k2*deltaT;
		k3=-wy*[Raux(:,3),knull,-Raux(:,1)];
		Raux=R+k3*deltaT;
		k4=-wy*[Raux(:,3),knull,-Raux(:,1)];
   FI=V1+(k1+2*k2+2*k3+k4)/6;
%%   fprintf('tiempo Vwake %f, tiempo Vbound %f, tiempo Vinf %f, tiempo Vrotante %f',t1,t2-t1,t3-t2,t4-t3)
elseif method==2   
 #Runge Kutta 2
		V1=AeroVinducida_rev5(R,Gw,wiB,wjB,wkB,wlB,1); #Velocidad inducida por estela
		V1=V1+AeroVinducida_rev5(R,G,niB,njB,nkB,nlB,1);	#Velocidad inducida por superficie
		V1=V1+repmat(Vlibre,nn,1)-wy*[R(:,3),knull,-R(:,1)]; #Velocidad libre y rotacion
    Raux=R+V1*deltaT;
    wiaux=Raux(conwake(:,1),:);
    wjaux=Raux(conwake(:,2),:);
    wkaux=Raux(conwake(:,3),:);
    wlaux=Raux(conwake(:,4),:);
		V2=AeroVinducida_rev5(Raux,Gw,wiaux,wjaux,wkaux,wlaux,1); #Velocidad inducida por estela
		V2=V2+AeroVinducida_rev5(Raux,G,niB,njB,nkB,nlB,1);	#Velocidad inducida por superficie
		V2=V2+repmat(Vlibre,nn,1)-wy*[Raux(:,3),knull,-Raux(:,1)]; #Velocidad libre y rotacion   
    
    FI=(V1+V2)/2;
elseif method==4    
    #Runge Kutta 4
		V1=AeroVinducida_rev5(R,Gw,wiB,wjB,wkB,wlB,1); #Velocidad inducida por estela
		V1=V1+AeroVinducida_rev5(R,G,niB,njB,nkB,nlB,1);	#Velocidad inducida por superficie
		V1=V1+repmat(Vlibre,nn,1)-wy*[R(:,3),knull,-R(:,1)]; #Velocidad libre y rotacion
%		
    Raux=R+0.5*V1*deltaT;
    wiaux=Raux(conwake(:,1),:);
    wjaux=Raux(conwake(:,2),:);
    wkaux=Raux(conwake(:,3),:);
    wlaux=Raux(conwake(:,4),:);
    
  	V2=AeroVinducida_rev5(Raux,Gw,wiaux,wjaux,wkaux,wlaux,1); #Velocidad inducida por estela
		V2=V2+AeroVinducida_rev5(Raux,G,niB,njB,nkB,nlB,1);	#Velocidad inducida por superficie
		V2=V2+repmat(Vlibre,nn,1)-wy*[Raux(:,3),knull,-Raux(:,1)]; #Velocidad libre y rotacion   
%    
    Raux=R+0.5*V2*deltaT;
    wiaux=Raux(conwake(:,1),:);
    wjaux=Raux(conwake(:,2),:);
    wkaux=Raux(conwake(:,3),:);
    wlaux=Raux(conwake(:,4),:);
%    
   	V3=AeroVinducida_rev5(Raux,Gw,wiaux,wjaux,wkaux,wlaux,1); #Velocidad inducida por estela
		V3=V3+AeroVinducida_rev5(Raux,G,niB,njB,nkB,nlB,1);	#Velocidad inducida por superficie
		V3=V3+repmat(Vlibre,nn,1)-wy*[Raux(:,3),knull,-Raux(:,1)]; #Velocidad libre y rotacion   
%    
    Raux=R+V3*deltaT;
    wiaux=Raux(conwake(:,1),:);
    wjaux=Raux(conwake(:,2),:);
    wkaux=Raux(conwake(:,3),:);
    wlaux=Raux(conwake(:,4),:);
    
	  V4=AeroVinducida_rev5(Raux,Gw,wiaux,wjaux,wkaux,wlaux,1); #Velocidad inducida por estela
		V4=V4+AeroVinducida_rev5(Raux,G,niB,njB,nkB,nlB,1);	#Velocidad inducida por superficie
		V4=V4+repmat(Vlibre,nn,1)-wy*[Raux(:,3),knull,-Raux(:,1)]; #Velocidad libre y rotacion   
  
    FI=(V1+2*V2+2*V3+V4)/6;

    endif

endfunction 
