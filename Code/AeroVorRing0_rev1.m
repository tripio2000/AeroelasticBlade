function [V]=AeroVorRing0_rev1(G,Ri,Rj,Rk,Rl,Rcp,La)
## Calcula velocidades inducidas por cada anillo vorticoso
## formados por los nodos Ri,Rj,Rk,Rl en los puntos Rcp 
## rev1: 5/12/18 incluyo tamaño del anillo "la" que es el promedio de las diagonales
	
#Ra=[Ri;Rj;Rk;Rl];#inicio segmento,tama�o (4mn,3)
#Rb=[Rj;Rk;Rl;Ri];#fin segmento, tama�o (4mn,3)
#Rc=[Rcp;Rcp;Rcp;Rcp];#punto de evaluacion, tama�o (4mn,3)
#G=[G;G;G;G]; #tama�o (4mn,1)
#[Vvec]=vorseg_vec (G,Ra,Rb,Rc);#tama�o (4mn,3)
#tam=size(Ri,1);
#V=Vvec(1:tam,:)+Vvec(1+tam:2*tam,:)+Vvec(1+2*tam:3*tam,:)+Vvec(1+3*tam:4*tam,:);
#V es la velocidad inducida por los cuatro segmentos del anillo en rcp, tama�o (mn,3)

gpumode=1;

if gpumode==1
  #GPU processing
  V=AeroVorSegGPU(G,Ri,Rj,Rcp);
  V=V+AeroVorSegGPU(G,Rj,Rk,Rcp);
  V=V+AeroVorSegGPU(G,Rk,Rl,Rcp);
  V=V+AeroVorSegGPU(G,Rl,Ri,Rcp);
  else
  V=AeroVorSeg_rev1(G,Ri,Rj,Rcp);
  V=V+AeroVorSeg_rev1(G,Rj,Rk,Rcp);
  V=V+AeroVorSeg_rev1(G,Rk,Rl,Rcp);
  V=V+AeroVorSeg_rev1(G,Rl,Ri,Rcp);
endif

endfunction
