% Esercizio 1.2. Si modifichi la precedente routine cosi` da studiare
% l'approssimazione, similmente all'esempio precedente, delle funzioni:
% abs(x-0.3), exp(x^2), exp(x), sin(x), sinc(x) = sin(x)/x.

clear all; esempio = 6;
d = [-1, 1]; %intervallo

switch esempio
    case 1
        ff = @(x) 1./(1+25*x.^2); % Funzione di Runge
    case 2
        ff = @(x) abs(x-0.3);     % Funz. continua con punto angoloso
    case 3
        ff = @(x) exp(x.^2);      % Funz. analitica
    case 4
        ff = @(x) exp(x);         % Funz. analitica
    case 5
        ff = @(x) sin(x);         % Funz. analitica
    case 6
        ff = @(x) sin(pi*x)./(pi*x);      % Funz. sinc
end

f1 = chebfun(ff, d, 'splitting', 'on');      % funzione (splitting on)
f2 = chebfun(ff, d);                         % funzione

fc_98  = chebfun(ff, 98);                          % interp. su 98 nodi di Cheb.
fc_100 = chebfun(ff, 100);                         % interp. su 100 nodi di Cheb.

err_98_split  = norm(f1 - fc_98, inf);             % errore ~ C * gamma^98  (splitting on)
err_98  = norm(f2 - fc_98, inf);                   % errore ~ C * gamma^98
err_100_split = norm(f1 - fc_100, inf);            % errore ~ C * gamma^100 (splitting on)
err_100 = norm(f2 - fc_100, inf);                  % errore ~ C * gamma^100

gamma_split = sqrt(err_100_split / err_98_split);  % approx. di gamma (splitting on)
gamma = sqrt(err_100 / err_98);                    % approx. di gamma
C_100_split = err_100_split / (gamma_split^100);   % approx. di C (splitting on)
C_100 = err_100 / (gamma^100);                     % approx. di C

nome_funz = func2str(ff);                           
disp(['Funzione in esame: ' strrep(nome_funz(5:end),'.','')]);    % stampa l'espressione della funzione

fprintf('\nErrore in norma inf. (con ''splitting on''):     ~ %f * %f ^n \n', C_100_split, gamma_split);
fprintf('Errore in norma inf. (senza ''splitting on''):   ~ %f * %f ^n \n\n', C_100, gamma);

[p_100, err_remez] = minimax(ff,100);               % confronto con la miglior approx.
fprintf('Err. calcolato con Remez / Err. con la stima sopra (senza splitting): %f.\n', err_remez/err_100);
fprintf('Err. calcolato con Remez / Err. con la stima sopra (con lo splitting): %f.\n\n', err_remez/err_100);

fprintf('Punti di Chebyshev usati per approssimare   %s   alla precisione di macchina: \n',strrep(nome_funz(5:end),'.',''));
fprintf('Con ''splitting on'': %d \n', length(f1));
fprintf('Senza ''splitting on'': %d \n\n', length(f2));
