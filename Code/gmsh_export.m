## Exporta la malla en un archivo .MSH para abrir con el GMSH
## Author: Nico Tripp

function gmsh_export (nodos,elem,nodoswake,conwake,G,dataout,m,n,filename,tstep)
# filename='malla2.msh';
fid = fopen (filename,'w');
fprintf(fid,'$MeshFormat \n');
fprintf(fid,'2.2 0 8 \n');
fprintf(fid,'$EndMeshFormat \n');


fprintf(fid,'$Nodes \n');
fprintf(fid,"%i \n",size(nodos,1)+size(nodoswake,1));#numero de nodos
nnodos=size(nodos,1);
idnodos=1:nnodos;
A=[idnodos',nodos(:,1:3)];#[id nodo,X, Y, Z]
idnodoswake=nnodos+1:nnodos+size(nodoswake,1);
B=[idnodoswake',nodoswake(:,1:3)];#[id nodo,X, Y, Z]
C=[A;B];
dlmwrite(fid,C,' ');
fprintf(fid,'$EndNodes \n');

fprintf(fid,'$Elements \n');
fprintf(fid,'%i \n',size(elem,1)+size(conwake,1));#numero de elementos
numelems=size(elem,1);
idelem=1:numelems;
A=[idelem',3*ones(numelems,1),1*ones(numelems,1),1*ones(numelems,1),elem(:,1:4)];#[elm-number elm-type number-of-tags < tag > ... node-number-list]
numwake=size(conwake,1);
idelemwake=numelems+1:numelems+numwake;
B=[idelemwake',3*ones(numwake,1),1*ones(numwake,1),2*ones(numwake,1),conwake(:,1:4)+nnodos*ones(numwake,1:4)];#[elm-number elm-type number-of-tags < tag > ... node-number-list]
C=[A;B];
dlmwrite(fid,C,' ');
fprintf(fid,'$EndElements \n');


fprintf(fid,'$ElementData \n');
fprintf(fid,'%i \n',1);			#number-of-string-tags
fprintf(fid,'Vortex strength \n');	#< "string-tag" >
fprintf(fid,'%i \n',1);			#number-of-real-tags
fprintf(fid,'%i \n',0);	#< real-tag > tiempo
fprintf(fid,'%i \n',3);			#number-of-integer-tags
fprintf(fid,'%i \n',0);	#time step
fprintf(fid,'%i \n',1);			#components (1=scalar)
fprintf(fid,'%i \n',numelems+numwake);	#elementos
datos=[idelem',G;idelemwake',conwake(:,5)]; #elm-number value
dlmwrite(fid,datos,' ');

%$NodeData
%1 #number-of-string-tags
%"A vector view" #< "string-tag" >
%1 #number-of-real-tags
%0.0 	#< real-tag > tiempo
%3  #number-of-integer-tags
%0 #time step
%3 #components (3=vector)
%6 #six associated nodal values
%1 0.0 0.0 0.0
%2 0.1 0.1 0.1
%3 0.2 0.2 0.2
%4 0.0 0.0 0.0
%5 0.2 0.2 0.2
%6 0.4 0.4 0.4
%$EndNodeData


fprintf(fid,'$EndElementData \n');



fclose (fid);




endfunction
