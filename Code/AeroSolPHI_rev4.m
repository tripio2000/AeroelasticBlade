function [G,Vindvecfull]=AeroSolPHI_rev4(nodos,elem,VlibreB,Riwake,Rjwake,Rkwake,Rlwake,Gwake,numelm,arranque)
## Vectorizado completo 7-5-2014
## modularizado evaluacion de vinducida 7-5-2014
## todo se resuelve en el sistema B
## rev4: uso voring0_pafun

	# nodos(:,x,y,z,u,v,w)
	# elem(:,nodoi,nodoj,nodok,nodol,rcpx,rcpy,rcpz,nx,ny,nz,ucpx,ucpy,ucpz)
	
	ri=nodos(elem(:,1),1:3);
	rj=nodos(elem(:,2),1:3);
	rk=nodos(elem(:,3),1:3);
	rl=nodos(elem(:,4),1:3);
  rcp=elem(:,5:7);
	sizeestela=size(Gwake,1);
	normales=elem(:,8:10)';
	VScpi=elem(:,11:13); 	#Velocidad absoluta por el movimiento de la superficie = Vrelativa + Vorigen + w x r
	G1=ones(numelm,1);

#	'Evaluacion matriz influencia'
%	Avec=voring0_pafun(G1vec,rivec,rjvec,rkvec,rlvec,rcpvec);
  Avec=AeroVinducida_rev5(rcp,G1,ri,rj,rk,rl,0);
%	disp('Vbound en cp con G unitario para matriz influencia A')

  nrx=rows(rcp);    
	Avecx=reshape(Avec(:,1:nrx),numelm,numelm).*(elem(:,8));
	Avecy=reshape(Avec(:,nrx+1:2*nrx),numelm,numelm).*(elem(:,9));
	Avecz=reshape(Avec(:,2*nrx+1:end),numelm,numelm).*(elem(:,10));
	A=Avecx+Avecy+Avecz;

  
	#Calculo el vector del lado derecho
	VlibreBvec=repmat(VlibreB,numelm,1);
	if arranque==0 #hay estela, incluyo calculo de velocidades
		[Vindvecfull]=AeroVinducida_rev5(elem(:,5:7),Gwake,Riwake,Rjwake,Rkwake,Rlwake,1);
%    disp('Vwake en cp para solphi')
		RHSx(:,1)=(VlibreBvec(:,1)+Vindvecfull(:,1)-VScpi(:,1)).*elem(:,8);
		RHSy(:,1)=(VlibreBvec(:,2)+Vindvecfull(:,2)-VScpi(:,2)).*elem(:,9);
		RHSz(:,1)=(VlibreBvec(:,3)+Vindvecfull(:,3)-VScpi(:,3)).*elem(:,10);
	else #arranque==1 no hay estela
    Vindvecfull=zeros(rows(elem),3);
		RHSx(:,1)=(VlibreBvec(:,1)-VScpi(:,1)).*elem(:,8);
		RHSy(:,1)=(VlibreBvec(:,2)-VScpi(:,2)).*elem(:,9);
		RHSz(:,1)=(VlibreBvec(:,3)-VScpi(:,3)).*elem(:,10);
	endif	
	RHSvec=-[RHSx+RHSy+RHSz];

	#Resuelvo el sistema
	G(:,1)=A\RHSvec;
  
%  size(A)
%  size(RHSvec)
%  size(G)
%  
endfunction	
