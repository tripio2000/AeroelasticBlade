## Exporta la malla en un archivo .MSH para abrir con el GMSH
## Author: Nico Tripp

function gmsh_export_bladecp(nodos,elem,G,dataout,m,n,filename,tstep)
# filename='malla2.msh';
fid = fopen (filename,'w');
fprintf(fid,'$MeshFormat \n');
fprintf(fid,'2.2 0 8 \n');
fprintf(fid,'$EndMeshFormat \n');


fprintf(fid,'$Nodes \n');
fprintf(fid,"%i \n",size(nodos,1));#numero de nodos
nnodos=size(nodos,1);
idnodos=1:nnodos;
A=[idnodos',nodos(:,1:3)];#[id nodo,X, Y, Z]
dlmwrite(fid,A,' ');
fprintf(fid,'$EndNodes \n');


fprintf(fid,'$Elements \n');
fprintf(fid,'%i \n',size(elem,1));#numero de elementos
numelems=size(elem,1);
idelem=1:numelems;
A=[idelem',3*ones(numelems,1),1*ones(numelems,1),1*ones(numelems,1),elem(:,1:4)];#[elm-number elm-type number-of-tags < tag > ... node-number-list]
dlmwrite(fid,A,' ');
fprintf(fid,'$EndElements \n');

fprintf(fid,'$ElementData \n');
fprintf(fid,'%i \n',1);			#number-of-string-tags
fprintf(fid,'Vortex strength \n');	#< "string-tag" >
fprintf(fid,'%i \n',1);			#number-of-real-tags
fprintf(fid,'%i \n',0);	#< real-tag > tiempo
fprintf(fid,'%i \n',3);			#number-of-integer-tags
fprintf(fid,'%i \n',0);	#time step
fprintf(fid,'%i \n',1);			#components (1=scalar)
fprintf(fid,'%i \n',numelems);	#elementos
cpelem=dataout(:,5); #Adimensional
datos=[idelem',cpelem]; #elm-number value
dlmwrite(fid,datos,' ');
fprintf(fid,'$EndElementData \n');



fclose (fid);




endfunction
