% Esercizio 1.3. Si modifichi la routine esempio1.m utilizzando
% quanto visto in esempio2.m, cosi` da paragonare la miglior approssimante
% di grado n = 10, 20, ..., 50 con l'interpolante in n+1 nodi di Chebychev,
%  su abs(x-0.3), exp(x^2), exp(x), sin(x), sinc(x) = sin(x)/x, cos(x),
%  tan(x)

clear all; esempio = 1;
d = [-1, 1]; % intervallo

switch esempio
    case 1
        ff = @(x) 1./(1+25*x.^2); % Funzione di Runge
    case 2
        ff = @(x) abs(x-0.3); % Funz. continua con punto angoloso
    case 3
        ff = @(x) exp(x.^2); % Funz. analitica
    case 4
        ff = @(x) exp(x); % Funz. analitica
    case 5
        ff = @(x) sin(x); % Funz. analitica
    case 6
        ff = @(x) sinc(x); % Seno normalizzato
    case 7
        ff = @(x) cos(x); % Funz. analitica
    case 8
        ff = @(x) tan(x); % Funz. analitica nell'intervallo dato
    case 9
        ff = @(x) atan(x); % Funz. analitica
    case 10
        ff = @(x) cosh(x); % Funz. analitica

end

nome_funz = func2str(ff);
disp(['Funzione in esame: ' strrep(nome_funz(5:end),'.','')]);    % stampa l'espressione della funzione

f1 = chebfun(ff,d,'splitting','on');    % chebfun della funzione (con 'splitting on')
f2 = chebfun(ff,d);                     % chebfun della funzione (senza 'splitting on')

L1 = length(f1); L2 = length(f2);
fprintf("\nPunti usati dall'interpolazione di Chebychev per raggiungere la precisione di macchina\n")
fprintf("- con splitting on: %f.\n- senza splittig on: %f\n\n",L1,L2)

if L1 == L2
    fprintf('Quindi non serve lo ''splitting on'', e ci possiamo fermare proprio a %d iterazioni.\n',L2);
    nn = 10:1:L2;

else
    fprintf('Quindi serve lo ''splitting on''.\n');
    if L2 > 200
        fprintf("Il numero di punti risulta %d. Fermiamoci a 200.\n",L2);
    end
    nn = 10 : 20 : 200;
    
end

fprintf('\nQui n = (#nodi - 1), o il grado della miglior approssimante.\n')

% inizializzazione vettori errori
err_cheb = [];        % Chebychev 
err_remez = [];       % Remez

% tabella errori
fprintf('\n n\t\t \t CHEBYCHEV \t\t     REMEZ \t\t \t punti usati Cheb. \t\t punti usati Remez \n');

for n=nn
    fc_n = chebfun(ff, n+1);      % interp. di Cheb. in n+1 nodi
    p_n = minimax(ff, n);         % miglior approx. di grado n
    
    err_cheb = [err_cheb norm(f1 - fc_n, inf)];
    err_remez = [err_remez norm(f1 - p_n, inf)];
    
    format long
    fprintf("%d\t \t \t %.5e \t\t %.5e  \t %d\t\t\t\t\t\t  %d\n",n,err_cheb(end),err_remez(end),length(fc_n),length(p_n));
end

semilogy(nn,err_cheb,'r--',nn,err_remez,'b--');
legend('Chebychev (rosso)','Remez (blu)');
title(strrep(nome_funz(5:end),'.',''))


