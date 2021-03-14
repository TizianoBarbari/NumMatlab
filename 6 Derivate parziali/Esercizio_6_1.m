clear all; hold off;
% Esercizio 6.1
% Si risolva l’equazione di Poisson nel quadrato [0,pi] x [0,pi] 
% con le date condizioni al contorno.

    f=inline('zeros(size(x))','x','y');
    g_down=inline('zeros(size(x))','x','y');
    g_up=inline('zeros(size(x))','x','y');
    g_left=inline('sin(y)','x','y');
    g_right=inline('(exp(pi)).*sin(y)','x','y');
    
    solution=inline('(exp(x)).*sin(y)','x','y');

nn=[2 4 8 16 32]
for n=1:length(nn)
    h=1/(nn(n)+1); x=pi*(h:h:1-h)'; y=x;
    [X,Y]=meshgrid(x,y);
    X=X'; Y=Y';

    u=poisson5pts(nn(n),f,g_left,g_right,g_down,g_up);
    
    % Uso reshape perché la soluzione abbia le stesse dimensioni di X e Y.
    Z=reshape(u,nn(n),nn(n));
    
   
        V=feval(solution,X,Y);
        
        err(nn(n))=norm(V(:)-Z(:),inf);
        if nn(n) == 1
            fprintf('\n \t [n]: %4.0f [ERR]: %2.2e',nn(n),err(nn(n)));
        else
            fprintf('\n \t [n]: %4.0f [ERR]: %2.2e [RATIO]: %2.2f',nn(n),err(nn(n)),err(nn(n)-1)/err(nn(n)));
        end
    
figure(n)
surf(X,Y,Z);                                % GRAFICO APPROSSIMATO
titolo = strcat('Metodo alle differenze con n = ',num2str(nn(n)));
title(titolo);
end
fprintf('\n');


hh=1/50; xx=pi.*(hh:hh:1-hh)'; yy=xx;
[XX,YY]=meshgrid(xx,yy);
XX=XX'; YY=YY';
figure(length(nn)+1)
surf(XX,YY,solution(XX,YY)); hold on;                 % GRAFICO SOLUZIONE ESATTA
title('Grafico della soluzione esatta');