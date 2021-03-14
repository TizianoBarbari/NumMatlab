clear all;

% Esercizio 4.1. Data la matrice di Hilbert di ordine 5,
% si calcolino col metodo delle potenze i suoi minimi e massimi autovalori
% in modulo.
% Da questi si determini il condizionamento della matrice in norma 2 e lo si confronti
% con cond(hilb(5),2). Eseguire lo stesso esercizio utilizzando il metodo QR.

A = hilb(5);             % matrice di Hilbert
toll = 10e-8;            % tolleranza
z0 = ones(size(A,1),1);  % vettore iniziale
nmax = 200;              % num. max. di iterazioni

fprintf('Si vogliono calcolare gli autovalori (di moduli minimo e massimo) della matrice di Hilbert di ordine 5.\n\n');
fprintf('- Metodo delle potenze -\n');

% AUTOVALORE DI MASSIMO MODULO
[lambda_max, v_max, n_iter_max, err_max]= metodo_potenze(A,z0,toll,nmax);
fprintf('\nDopo %d interazioni:\n',n_iter_max);
fprintf('Autovalore di mod. max.: %f\n',lambda_max(end));

% AUTOVALORE DI MINIMO MODULO
[lambda_min, v_min, n_iter_min, err_min]= metodo_potenze(inv(A),z0,toll,nmax);
fprintf('\nDopo %d interazioni:\n',n_iter_min);
fprintf('Autovalore di mod. min.: %f\n',1/lambda_min(end));
% Il metodo delle potenze applicato a A^(-1) restituisce il reciproco
% dell'autovalore di modulo minimo di A.

cond1=cond(hilb(5),2); cond2=lambda_max(end)*lambda_min(end);
fprintf('\nCondizionamento (norma eucl.) di Hilb(5) calcolato in due modi:\n');
fprintf('1) Routine di MATLAB: cond(hilb(5),2) =\t %7.10f.\n',cond1);
fprintf('2) max |lambda_i| / min |lambda_i| =\t %7.10f.\n',cond2);
fprintf('Discrepanza: %7.10f = %1.10e.\n\n',cond1-cond2,cond1-cond2);

% ----------------------------------------------------------
% FATTORIZZAZIONE QR

[H,Q]=houshess(A);
format long
[T,hist]=metodo_QR(H, nmax);
fprintf(' - Metodo QR - \n\n');
t= diag(T); disp(t); cond3=t(1)/t(end);
fprintf('\nDunque:\nlambda_max = %7.10f\nlambda_min = %7.10f.\n',t(1),t(end));
fprintf('\nCon tali valori cond(hilb(5),2) = %7.10f \n',cond3);
fprintf('mentre la routine di MATLAB porge %7.10f \n',cond1);
fprintf('Discrepanza: %7.10f = %1.10e\n',cond1-cond3,cond1-cond3);



