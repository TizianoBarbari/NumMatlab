% Esercizio 1.4. Verificare che per nodi equispaziati, unisolventi a grado n,
% per n = 1, ..., 50 si ha che 
% (2^(n-2))/(n^2)   <=  LebesgueConstant   <=   (2^(n+3))/n,
% come stabilito in Two results on Polynomial Interpolation in Equally Spaced Points,
% da Trefethen e Weideman, J.A.T. (65), 247-260 (1991).

% E` buona la stima asintotica Lamda_n ~ (2^(n+1))/(exp(1)*n*ln(n)) ?

clear all;
fprintf('\n \t COSTANTI DI LEBESGUE, NODI EQUISPAZIATI \n');
max_deg=50;
costanti_di_lebesgue=zeros(max_deg,1);
stima_1=zeros(max_deg,1);
stima_2=zeros(max_deg,1);
stima_asint=zeros(max_deg,1);
fprintf('\nGRADO \t STIMA DAL BASSO \t COST. DI LEB.\t STIMA DALL''ALTO \n');

nn=1:max_deg;
for deg=nn
    x=linspace(-1,1,deg+1);
    [L,L_const]=lebesgue(x); % FUNZIONE E COSTANTE DI LEBESGUE (MATLAB routine)
    lebesgue_constants(deg) = L_const;
    stima_1(deg) =  (2^(deg-2))/(deg^2) ; % stima dal basso cost.leb.
    stima_2(deg)  =  (2^(deg+3))/deg;      % stima dall'alto cost.leb.
    fprintf('%3.0f         %1.4f         %1.4e         %1.4e \n',deg,stima_1(deg),L_const,stima_2(deg));
    
    stima_asint(deg) = (2^(deg+1))/(exp(1)*deg*log(deg));
end

semilogy(nn,stima_1,'b-',nn,lebesgue_constants,'r-',nn,stima_2,'p-',nn,stima_asint,'ro');
title('Stima della costante di Lebesgue per nodi equispaziati');
[aa,bb]=legend('$\frac{{2}^{n-2}}{n ^2}$ (Stima dal basso)','$\Lambda_n$ (MATLAB routine)','$\frac{{2}^{n+3}}{n}$ (Stima dall''alto)','Stima $\Lambda_n \sim \frac{{2}^{n+1}}{exp(1)*n*ln(n)}$','Interpreter','latex','location','northwest');