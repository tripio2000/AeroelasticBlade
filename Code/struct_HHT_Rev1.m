function [Znuevo,ddqdtnuevo]=struct_HHT_Rev1(qactual,dqdtactual,M,K,fnuevo,factual,DT,TC)
%Solucion por metodo implicito HHT, alfa atenúa la respuesta a frecuencias superiores a 1/(2h)
%qactual: coordenadas generalizadas actuales
%dqdtactual: velocidades generalizadas actuales
%factual: fuerzas generalizadas actuales
%fnuevo: fuerzas generalizadas siguientes
%Las ecuaciones están en formato adimensional (en el tiempo) por TC
###Rev1 9/3/19: corrijo calculo de aceleracion actual

alfa=0.02; %entre 0 a 1/3 elegido 0.02
bet=((1+alfa)^2)/4; %parámetro para posicion
gama=1/2+alfa; %parámetro para velocidad
  
%Estado actual
ddqdtactual=M\((TC^2)*(factual-K*qactual)); %aceleracion actual
ddqdtnuevo=(M+(TC^2)*K*DT^2*bet*(1-alfa))\((TC^2)*((1-alfa)*fnuevo+alfa*factual-K*(qactual+(1-alfa)*DT*dqdtactual+(1-alfa)*DT^2*(1/2-bet)*ddqdtactual))); %prediccion aceleracion
dqdtnuevo=dqdtactual+DT*((1-gama)*ddqdtactual+gama*ddqdtnuevo); %prediccion velocidad
qnuevo=qactual+DT*dqdtactual+DT^2*((1/2-bet)*ddqdtactual+bet*ddqdtnuevo); %prediccion posicion
  
Znuevo=[qnuevo;dqdtnuevo]; %salida de datos
  
endfunction
