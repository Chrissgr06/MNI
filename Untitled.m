% f = @(c) (667.38./c)*(1-exp(-0.146843.*c))-40;
% %Condiciones iniciales
% tol = 0.01;
% ea=100;
% 
% xl=12;
% xu=16;
% %xra = 16;
% 
% %Procedimiento
% while ea > tol
%     xr = (xl+xu)/2;
%    % ea = abs((xr - xra)/xr)*100
%     if f(xl)*f(xr) < 0
%         xu = xr;
%     end
%     if f(xr)*f(xu) < 0
%         xl = xr;
%     end
%     xra = xr
% end
%-----------------------------------
%Punto Fijo
% xo = 10;
% tol = 0.01;
% imax = 5;
% 
% g1 = @(x) x - (x^3) - 4*(x^2) + 10;
% g2 = @(x)((10/x)-(4*x))^(1/2);
% g3 = @(x)1/2*(10 - (x^3))^(1/2);
% g4 = @(x)(10/(4 + x))^(1/2);
% g5 = @(x)x - (((x^3) + (4*(x^2)) - 10)/((3*x^2)+8));
% 
% for i=1:imax
%     x = g(xo)
%     eo = (x-xo)/x
%     if ea <tol
%         x
%         ea
%     else
%         'no hay convergencia'
%     end
%     
% end
%------------------------------------
% %N-R Method
% % h1 = @(x) ((x^3) + (4*(x^2)) - 10); xo = 1;
% % dh1 = @(x) (3*(x^2) + (8*x));
% % % h2 = @(x) (exp(-x)) - x; xo = 0;
% % % h3 = @(x) (x^10)-1; xo = 0.5;
% % 
% % tol = 0.1;
% % imax = 100;
% % tic;
% % for i=1:imax
% %     x = xo -(h1(xo)/dh1(xo))
% %     ea = abs((x-xo)/x)
% %     if ea<tol
% %         break
% %     else
% %         xo=x;
% %     end
% % end
% % toc;
% 
% %--------------------------
% %Regla de cramer
% 
% %Matriz de coordenadas generalizadas
% a = [0.3,0.52,1;0.5,1,1.9;0.1,0.3,0.5];
% %Vector de terminos independientes
% b = [-0.01,0.67,-0.44];
% %Determinantes
% d = det(a);
% %Solución
% a1=a;
% a1(:,1) = b;
% x1=det(a1)/d
% 
% a2=a;
% a2(:,2) = b;
% x2=det(a2)/d
% 
% a3=a;
% a3(:,3) = b;
% x3=det(a3)/d
%-----------------------------------
%Método de Gauss
%Matriz de coordenadas generalizadas
a = [0.3,0.52,1;0.5,1,1.9;0.1,0.3,0.5];
%Vector de terminos independientes
b = [-0.01,0.67,-0.44];
function [a,b] = Eliminacion(a,b)
n=size(a,1);
%Ciclo de pivote
for p=1:n-1
    %Ciclo de filas
    for i=p+1:n
        %Factor de pivotaje: elemento pivote = a(p,p)
        factor = a(i,p)/a(p,p);
        %ciclo de columnas (incognitas)
        for j=p:n
            a(i,j) = a(i,j) - factor*a(p,j);
        end
        b(i) = b(i) - factor*b(p);
    end
end
end

function[x] = Sust(a,b)
n = size(a,1);
%Primer despeje
x(n) = b(n)/u(n,n);

%Ciclo por filas hacia atras
for i=n-a:-1:1
    suma=b(i);
    %Ciclopor columnas
    for j=i+1:n
        suma=suma-a(i,j)*x(j);
    end
    %Calcular x
    x(i) = suma/a(i,i)
end
end




