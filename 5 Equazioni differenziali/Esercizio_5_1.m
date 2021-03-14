clear all; hold off;
% Esercizio 5.1. Si risolva utilizzando il metodo di Eulero esplicito, con passi h uguali
% rispettivamente a 0.2, 0.1, 0.05, 
%       y'(x) = (cos(y(x)))^2    [0 <= x <= 10]
%       y(0) = 0,
% la cui soluzione e` y(x) = arctan(x).

% -------------------------------------------
% DATI DEL PROBLEMA DI CAUCHY
y0=0; t0=0; t_end=10; passi=[0.2 0.1 0.05];
funz = inline('cos(x).^2');

L=length(passi); M=(t_end-t0)/min(passi)+1;     % dimensioni matrice y (ogni colonna per un passo h)
y = zeros(M,L); x = zeros(M,L);                 % inizializza a 0 le matrici
y(1,:) = y0;                                    % y(0)=y0 (prima riga della matrice)

for h=1:L                                            % per ogni passo h si esegue il metodo di Eulero Esplicito
    n(h) = ((t_end-t0)/passi(h)+1)';                 % numero di punti
    x(1:n(h),h) = t0 : passi(h) : t_end;             % punti
    [t,yy]= eulero_esplicito(t0,y0,t_end,passi(h),funz);
    y(1:n(h),h) = yy;
end

figure(1)
plot(x(1:n(1),1),y(1:n(1),1),'ro'); hold on;
plot(x(1:n(2),2),y(1:n(2),2),'gp');
plot(x(1:n(3),3),y(1:n(3),3),'k*');
plot(x(:,3),atan(x(:,3)),'b-');    % soluzione vera
title('Soluzione numerica del problema di Cauchy y''(x) = (cos(y(x)))^2 [0 <= x <= 10], y(0) = 0');
legend('Eul. Espl. passo 0.2','Eul. Espl. passo 0.1','Eul. Espl. passo 0.05','Soluzione esatta (arctan(x))','location','southeast');
hold off;

% -------------------------------------------------
% ERRORI (ASSOLUTO e RELATIVO)
err_ass = zeros(M,L); err_rel = zeros(M,L);
for h=1:L
err_ass(1:n(h),h) = abs(atan(x(1:n(h),h))-y(1:n(h),h));
err_rel(1:n(h),h)= err_ass(1:n(h),h) ./ atan(x(1:n(h),h));
end

figure(2)
semilogy(x(1:n(1),1),err_ass(1:n(1),1),'ro',x(1:n(1),1),err_rel(1:n(1),1),'r*'); hold on;
semilogy(x(1:n(2),2),err_ass(1:n(2),2),'go',x(1:n(2),2),err_rel(1:n(2),2),'g*');
semilogy(x(1:n(3),3),err_ass(1:n(3),3),'ko',x(1:n(3),3),err_rel(1:n(3),3),'k*');
title('Errore di Eulero Esplicito');
legend('Err. ass. passo 0.2','Err. rel. passo 0.2','Err. ass. passo 0.1','Err. rel. passo 0.1','Err. ass. passo 0.05','Err. rel. passo 0.05');
hold off;

% -----------------------------------------------------------
% Per selezionati valori di x, ad esempio x = 10,
% si calcoli il rapporto con cui l’errore decresce quando h e`
% dimezzato.
% 

hh=[0.4]; for j=2:10 hh(j)=hh(j-1)/2; end       % 10 passi h dimezzati progressivamente
errass = zeros(length(hh),1);
for j=1:length(hh)
[t,yy] = eulero_esplicito(t0,y0,t_end,hh(j),funz);
errass(j) = abs(yy(end)-atan(10));
end
rapporti=errass(2:end)./errass(1:end-1);

figure(3);
plot(1:length(hh)-1,rapporti);
title('Rapporto con cui decresce l''errore per x=10');

fprintf('Per x=10, il rapporto con cui descresce l''errore, risulta: \n');
for j=1:length(hh)-1
fprintf('Da passo %f a passo %f:\t %f \n',hh(j), hh(j+1),rapporti(j));
end