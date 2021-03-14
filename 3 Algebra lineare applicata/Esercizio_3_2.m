% Si calcoli la matrice di Poisson P_20 di ordine 20.
% Sia b il vettore composto di componenti uguali a 1, avente lo stesso numero di righe di P_20.
% Si risolva col metodo del gradiente coniugato il problema P_20 x = b,
% con tolleranza di 10^(-12), partendo da x_0 = [0 ... 0]. Converge?

clear all; max_it=200; tol=10^(-12); siz=20;
A = makefish(siz);     % MATRICE DI POISSON.
b=ones(size(A,1),1);   % TERMINE NOTO.

x_sol=A\b;             % SOLUZIONE ESATTA (METODO LU)

norm_x_sol=norm(x_sol);
if norm(x_sol) == 0
    norm_x_sol=1;
end

x=zeros(size(b));      % VALORE INIZIALE.

fprintf('MATRICE DI POISSON DI ORDINE 20: \n');

% GRADIENTE CONIUGATO.
M=eye(size(A));
[x_gc, error_gc, iter_gc, flag_gc] = cg(A, x, b, M, max_it, tol);

fprintf('\t \n [GRA.CON.] [STEP REL., NORMA 2]: %2.2e [REL.ERR.]: %2.2e',error_gc,norm(x_gc-x_sol)/norm_x_sol);
fprintf('\t \n            [ITER.]: %3.0f [FLAG]: %1.0f \n',iter_gc,flag_gc);

% JACOBI.
[x_j, error_j, iter_j, flag_j]  = jacobi(A, x, b, max_it, tol);

fprintf('\t \n [JACOBI  ] [STEP REL., NORMA 2]: %2.2e [REL.ERR.]: %2.2e',error_j,norm(x_j-x_sol)/norm_x_sol);
fprintf('\t \n            [ITER.]: %3.0f [FLAG]: %1.0f \n',iter_j,flag_j);

% GAUSS-SEIDEL.
w=1;
[x_gs, error_gs, iter_gs, flag_gs]  = sor(A, x, b, w, max_it, tol);

fprintf('\t \n [GAU.SEI.] [STEP REL., NORMA 2]: %2.2e [REL.ERR.]: %2.2e',error_gs,norm(x_gs-x_sol)/norm_x_sol);
fprintf('\t \n            [ITER.]: %3.0f [FLAG]: %1.0f \n',iter_gs,flag_gs);

%-----------------------------------------------------------------------
% Applicare il gradiente coniugato alla matrice minij di ordine 100
% e paragonare i risultati con quelli ottenuti da Jacobi, Gauss-Seidel e SOR ottimale.

clear all;
A = gallery('minij',100);
max_it=200; tol=10^(-12); siz=20;
b=ones(size(A,1),1);   % TERMINE NOTO.
x_sol=A\b;             % SOLUZIONE ESATTA (METODO LU)
x=zeros(size(b));      % VALORE INIZIALE.

norm_x_sol=norm(x_sol);
if norm(x_sol) == 0
    norm_x_sol=1;
end

fprintf('\n MATRICE minij DI ORDINE 100: \n');

% GRADIENTE CONIUGATO.
M=eye(size(A));
[x_gc, error_gc, iter_gc, flag_gc] = cg(A, x, b, M, max_it, tol);

fprintf('\t \n [GRA.CON.] [STEP REL., NORMA 2]: %2.2e [REL.ERR.]: %2.2e',error_gc,norm(x_gc-x_sol)/norm_x_sol);
fprintf('\t \n            [ITER.]: %3.0f [FLAG]: %1.0f \n',iter_gc,flag_gc);

% JACOBI.
[x_j, error_j, iter_j, flag_j]  = jacobi(A, x, b, max_it, tol);

fprintf('\t \n [JACOBI  ] [STEP REL., NORMA 2]: %2.2e [REL.ERR.]: %2.2e',error_j,norm(x_j-x_sol)/norm_x_sol);
fprintf('\t \n            [ITER.]: %3.0f [FLAG]: %1.0f \n',iter_j,flag_j);

% GAUSS-SEIDEL.
w=1;
[x_gs, error_gs, iter_gs, flag_gs]  = sor(A, x, b, w, max_it, tol);

fprintf('\t \n [GAU.SEI.] [STEP REL., NORMA 2]: %2.2e [REL.ERR.]: %2.2e',error_gs,norm(x_gs-x_sol)/norm_x_sol);
fprintf('\t \n            [ITER.]: %3.0f [FLAG]: %1.0f \n',iter_gs,flag_gs);

% SOR.
w_vett=0.8:0.025:2;

for index=1:length(w_vett)
    w=w_vett(index);
   [x_sor, error_sor(index), iter_sor(index), flag_sor(index)]  = sor(A, x, b, w, max_it, tol);
   relerr(index)=norm(x_sor-x_sol)/norm_x_sol;
end

[min_iter_sor, min_index]=min(iter_sor);

fprintf('\t \n [SOR OTT.] [STEP REL., NORMA 2]: %2.2e [REL.ERR.]: %2.2e',error_sor(min_index),relerr(min_index));
fprintf('\t \n            [ITER.]: %3.0f [FLAG]: %1.0f [w]: %2.3f \n',min_iter_sor,flag_sor(min_index),w_vett(min_index));
