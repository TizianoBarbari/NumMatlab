% Esercizio 1.1. Sapendo che l'errore in norma infinito dell'interpolante
% nei nodi di Chebyshev e` asintoticamente del tipo C*(gamma^n), determinare
% C e gamma.

clear all;
ff = @(x) 1./(1+25*x.^2);             % Funzione di Runge
f=chebfun(ff, [-1, 1], 'splitting', 'on'); 

fc_98  = chebfun(ff, 98);             % interp. su 98 nodi di Cheb.
fc_100 = chebfun(ff, 100);            % interp. su 100 nodi di Cheb.

err_98  = norm(f - fc_98, inf);       % errore ~ C * gamma^98
err_100 = norm(f - fc_100, inf);      % errore ~ C * gamma^100

gamma = sqrt(err_100 / err_98);       % approx. di gamma
C = err_100 / (gamma^100);            % approx. di C, usando quella di gamma

fprintf('Errore in norma infinito ~ C*gamma^n:\n\t %f * %f^n \n\n', C, gamma);
fprintf('Stima di gamma:\n\t %f.\n\n',gamma);
fprintf('Stima di C:\n\t %f.\n\n',C);

[p_99, err_remez] = minimax(f,99);  % confronto con miglior approx.
format long;
fprintf('Err. calcolato con Remez / Err. con la stima C*(gamma^n):\n\t %f.\n\n', err_remez/err_100);