% Usando una opportuna formula gaussiana, con n = 10, 20, 30 punti,
% calcolare integrali di diverse funzioni.
% Descrivere come decresce l’errore in scala semilogaritmica.

clear all; scelta = 4;

switch scelta
    case 1
        ff='exp(x).*cos(4*x)'; f = inline(ff); I = [0, pi];
        alpha = 0; beta = 0; true_val = (exp(pi)-1)/17;
    case 2
        ff='x.^(5/2)'; f=inline(ff); I = [0, 1];
        alpha = 0; beta = 0; true_val = 2/7;
    case 3
        ff='2^(-7/2)*(1+x).^0'; f = inline(ff); I = [-1, 1];
        % stessa funzione del caso 2 ma riscalata e riscritta
        % per usare il peso di Jacobi
        alpha = 0.; beta = 2.5; true_val = 2/7;
    case 4
        ff='exp(cos(x))'; f = inline(ff); I = [-pi, pi];
        alpha = 0; beta = 0; true_val = 7.95492652101284;
end

nome_funz = strrep(ff,'.','');

errori = [];
fprintf('Formula di Gauss-Jacobi con 10,20,30 punti.\n');
fprintf('\nFunzione da integrare: %s, tra %1.1f e %1.1f, usando il peso (1-x)^%1.1f (1+x)^%1.1f \n',nome_funz,I(1),I(2),alpha,beta)
fprintf('Risultato ''esatto'': %2.14f.\n\n',true_val);

pts=10:10:30;           % numero di punti
for p=pts
[I_val,x_jac,w_jac]=integrazione_gauss_jacobi(p, alpha, beta, I(1), I(2), f);
err=abs(true_val - I_val);
fprintf('Valore calcolato con %d punti: %2.16f.\tErrore con %d punti: %1.5e\n', p,I_val,p, err);
errori=[errori err];
end


semilogy(pts, errori,'r-');

title(sprintf('Errore di Gauss-Jacobi su %s',nome_funz),'Interpreter', 'none');

