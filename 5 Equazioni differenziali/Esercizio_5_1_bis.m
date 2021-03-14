clear all;
% Esercizio 5.1. Si risolva utilizzando il metodo di Eulero esplicito, con passi h uguali
% rispettivamente a 0.2, 0.1, 0.05, 
%       y'(x) = (cos(y(x)))^2    [0 <= x <= 10]
%       y(0) = 0,
% la cui soluzione e` y(x) = arctan(x).

% -------------------------------------------
% DATI DEL PROBLEMA DI CAUCHY
t0=0; t_end=10; y0=0; passi = [0.2, 0.1, 0.05]; % passi


L=length(passi); M=(t_end-t0)/min(passi)+1;     % dimensioni matrice y (ogni colonna per un passo h)
y = zeros(M,L); x = zeros(M,L);                 % inizializza a 0 le matrici
y(1,:) = y0;                                    % y(0)=y0 (prima riga della matrice)

for h=1:L                                 % per ogni passo h si esegue il metodo di Eulero Esplicito
n(h) = ((t_end-t0)/passi(h)+1)';
x(1:n(h),h) = t0: passi(h) : t_end;  % punti

    for i=1:n(h)-1                        % 
        f = cos(y(i,h))^2;                % funzione y'
        y(i+1,h) = y(i,h) + passi(h)*f;   % applicazione del metodo vera e propria
    end                                   % 

end                                  % fine ciclo for relativo al passo h-esimo

figure(1)
plot(x(1:n(1),1),y(1:n(1),1),'ro'); hold on;
plot(x(1:n(2),2),y(1:n(2),2),'gp');
plot(x(1:n(3),3),y(1:n(3),3),'k*');
plot(x(:,3),atan(x(:,3)),'b-');    % soluzione vera
title('Soluzione del P.C. y''(x) = (cos(y(x)))^2 [0 <= x <= 10], y(0) = 0');
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

% -----------------------------------------------------------
% Per selezionati valori di x, ad esempio x = 10,
% si calcoli il rapporto con cui l’errore decresce quando h e`
% dimezzato.
% 
% fprintf('Per x=10, al dimezzare di h, l''errore descresce con il rapporto %f.\n',err_ass(n(1),1) / err_ass(n(2),2));
% err_ass(n(1),1) / err_ass(n(2),2)
% err_ass(n(2),2) / err_ass(n(3),3)