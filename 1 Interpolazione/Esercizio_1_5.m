% Esercizio 1.5 (Facoltativo). Che errore si compie in norma infinito e in norma 2,
% approssimando la funzione f(t) = log(0.001 + t), t nell'intervallo [0, 1], con
% 1) polinomi algebrici di grado 11, 31, 51 via Chebfun?
% 2) polinomi trigonometrici di grado 11, 31, 51 via Chebfun (nella
% versione col ’trunc’)?

clear all;
ff = @(x) (log(0.001)+x);                      % funzione
d=[0 1];                                       % intervallo
f = chebfun(ff,d);

err_2_pol = []; err_inf_pol = [];              % inizializzazione vettori errori
err_2_trig = []; err_inf_trig = [];

nn=12:20:52;                                   % num. punti = grado polinomi + 1
for n=nn
    pol =  chebfun(ff,d,n);                      % pol. algebrici
    trig = chebfun(ff,d,'trunc',n, 'trig');      % pol. trig.
    err_2_pol = [err_2_pol norm(f-pol,2)];
    err_2_trig = [err_2_trig norm(f-trig,2)];
    err_inf_pol = [err_inf_pol norm(f-pol,inf)];
    err_inf_trig = [err_inf_trig norm(f-trig,inf)];    
end

figure(1)
semilogy(nn,err_2_pol,'b-',nn,err_2_trig,'b*',nn,err_inf_pol,'r-',nn,err_inf_trig ,'r*');
title('Errore, in norma 2 e inf., interpolando di log(0.001+t) con polinomi algebrici e trigonometrici');
[aa,bb]=legend('Alg. norma 2','Trig. norma 2','Alg. norma inf.','Trig. norma inf.','location','east');

figure(2)
hold on;
plot(chebfun(ff,d,'trunc',12, 'trig'),'color','red')
plot(chebfun(ff,d,'trunc',32, 'trig'),'color','yellow')
plot(chebfun(ff,d,'trunc',52, 'trig'),'color','blue')
plot(chebfun(ff,d),'color','black')
legend('trig. grado 11 (rosso)','trig. grado 31 (giallo)','trig. grado 51 (blu)','Funzione ''vera'' (nero)','location','west')
title('Interpolazione trigonometrica (con ''trunc.'')')
hold off;

