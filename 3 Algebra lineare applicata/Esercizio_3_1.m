clear all; hold off;

% Esercizio 3.1. Si calcoli la matrice simmetrica e definita positiva
% min_ij(20) di ordine 20.

A = gallery('minij',20); % Matrice (simm. e def. pos.)
b=ones(size(A,1),1);     % Termine noto
x=zeros(size(b));        % Vettore iniziale

x_sol=A\b;               % Soluzione esatta (metodo LU)

norm_x_sol=norm(x_sol);
if norm(x_sol) == 0
    norm_x_sol=1;
end

% ------------------------------------------------------
% Si risolva col metodo di Jacobi il problema minij(20) x = b,
% con tolleranza di 10^(-6), partendo da x_0 = [0 ... 0]. Converge?

max_it=1000; tol=10^(-6);
[x_j, error_j, iter_j, flag_j]  = jacobi(A, x, b, max_it, tol);

fprintf('\t \n [JACOBI  ] [STEP REL., NORMA 2]: %2.2e [REL.ERR.]: %2.2e',error_j,norm(x_j-x_sol)/norm_x_sol);
fprintf('\t \n            [ITER.]: %3.0f [FLAG]: %1.0f \n',iter_j,flag_j);

% ------------------------------------------------------
% Si risolva col metodo di SOR con omega = 0.01 : 0.01 : 1.99
% il problema minij(20) x = b, con tolleranza di 10^(-6), partendo da
% x_0 = [0 ... 0].

w_vett = 0.01 : 0.01 : 1.99;

for index=1:length(w_vett)
    w=w_vett(index);
   [x_sor, error_sor(index), iter_sor(index), flag_sor(index)]  = sor(A, x, b, w, max_it, tol);
   relerr(index)=norm(x_sor-x_sol)/norm_x_sol;
end

% Converge?
% Eseguire il plot in scala semilogaritmica, avendo in ascisse omega
% e in ordinate il numero di iterazioni eseguite.
% Quale sembra un buon parametro da utilizzare?
% Calcolare il parametro omega che minimizza il numero di iterazioni svolte da
% SOR.

[min_iter_sor, min_index]=min(iter_sor);

fprintf('\t \n [SOR OTT.] [STEP REL., NORMA 2]: %2.2e [REL.ERR.]: %2.2e',error_sor(min_index),relerr(min_index));
fprintf('\t \n            [ITER.]: %3.0f [FLAG]: %1.0f [w]: %2.3f \n',min_iter_sor,flag_sor(min_index),w_vett(min_index));
 
plot(w_vett,iter_sor,'r-');
title('Numero di iterazioni in funzione del parametro w (metodo SOR)');
xlabel('w');
ylabel('Iterazioni')

