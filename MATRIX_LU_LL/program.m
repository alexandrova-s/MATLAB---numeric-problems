clear all;
clc;
close all; 
 
N=[3,10,20];
a=sqrt(2)/2;
x=sqrt(a*a + 1/2)-1;
rozmiar = 3;
 
%Zadanie 1
A = macierz(N,1);

%Zadanie 2
wykresDet(N, a, rozmiar);
figure();
wykresCond(N, a, rozmiar);

%zadanie 3
for i = 1:rozmiar
    M = macierz(N(i), 1);
    inv(M);
    LL_(M);
    LU_(M);
end

%zadanie 4 i 5
wskaznikiBl(N, rozmiar)


%Zadanie1 - Generacja Macierzy 
function A = macierz(N,x)
A = ones(N);

for i=1:N
    A(i:end,i)=A(i,i)*i*(4/9);
    A(i,i:end)=A(i,i);
    A(i:end,1)=x*(-2/3);
    A(1,i:end)=x*(-2/3);
    A(1,1)=x^2;
end

end

%Zadanie 2
%wyznaczenie parametrów
function P = parametryMacierzy(N,a)
amin = a-0.01;
amax = a+0.01;
j = 0;
P = ones(N);
rozmiar=1000;
a = linspace(amin, amax, rozmiar);

    for i=1:rozmiar
        x = sqrt(a(i)^2 + 1/2)-1; 
        A = macierz(N,x);
        
    j = j+1;
    P(1,j)=a(i);
    P(2,j)=det(A);
    P(3,j)=cond(A);
  
    end
    
end

%Wykres wyznacznika macierzy
function D = wykresDet(N,a,rozmiar)
 
    for i =1:rozmiar  
        P = parametryMacierzy(N(i),a);
        q=P(1,:);
        w=P(2,:);
        semilogy(q,w);
        hold on;
    end
    grid on;
    
    xlabel('Parametr alfa');
    ylabel('Wyznacznik macierzy');
    title('Wykres det od alfa');
    
    hold off;
end

%Wykres wzkaŸnika uwarunkowania
function C = wykresCond(N,a,rozmiar)
 
    for i =1:rozmiar
        P = parametryMacierzy(N(i),a);
        q=P(1,:);
        e=P(3,:);
        semilogy(q,e);
        hold on ;
    end
    grid on;
    
    xlabel('Parametr alfa');
    ylabel('WskaŸnik uwarunkowania macierzy');
    title('Wykres cond od alfa');
    
    hold off;
end

%Zadanie3
%Rozk³ad LU
function LU = LU_(A)
n = size(A);
[L, U, P] = lu(A);
I=eye(n);
K=ones(n);
B=ones(n);
for k = 1:n
    for i = 1:n 
        sumL=0;        
        for j=1:i-1        
            if(i>1)
            sumL=sumL+L(i,j)*B(j,k);
            end             
        end        
        B(i,k)=(I(i,k)-sumL)/L(i,i);        
    end    
end

for a = 1:n
    for b = n:-1:1     
        sumU=0;        
        for c = n:-1:b+1
           sumU=sumU+U(b,c)*K(c,a);
        end        
        K(b,a)=(B(b,a)- sumU)/U(b,b);       
    end    
end

LU = K*P;

end

%rozk³ad LL
function LL = LL_(A)
n = size(A);
L=chol(A);
Lt=L';
I=eye(n);
A=ones(n);
B=ones(n);
for k = 1:n
    for i = 1:n 
        sumLt=0;        
        for j=1:i-1         
            if(i>1)
            sumLt=sumLt+Lt(i,j)*B(j,k);
            end             
        end        
        B(i,k)=(I(i,k)-sumLt)/Lt(i,i);        
    end    
end

for a = 1:n
    for b = n:-1:1       
        sumL=0;        
        for c = n:-1:b+1
           sumL=sumL+L(b,c)*A(c,a);
        end        
        A(b,a)=(B(b,a)- sumL)/L(b,b);       
    end    
end
LL = A;
end

%Zadnaie 4 i 5
%Norma druga
function delta = norm2(A)

    delta = max(sqrt(eigs(A' * A)));

end

%Norma nieskoñczona
function delta = normInf(A)

    delta = max(sum(abs(A), 2));

end

%Wykresy wzkaŸników
function Z5 = wskaznikiBl(N, rozmiar)

    K = 22;
    
    zakresK = 0 : (K - 1);
    zakresX = 2.^(zakresK) / 300;
    LLdelta2 = zeros(1, K);
    LLdeltaInf = zeros(1, K);
    LUdelta2 = zeros(1, K);
    LUdeltaInf = zeros(1, K);
    invDelta2 = zeros(1, K);
    invdeltaInf = zeros(1, K);

    for i = 1:rozmiar  
        
      I = eye(N(i));

      for k = zakresK
        x = 2 ^ k / 300;
        A = macierz(N(i), x);

%Dla metody LL

        invA = LL_(A);
        M = A * invA - I;
        delta2 = norm2(M);
        deltaInf = normInf(M);
        delta2norm = norm(M, 2);
        deltaInfnorm = norm(M, inf);
        LLdelta2(k + 1) = delta2;
        LLdeltaInf(k + 1) = deltaInf;
        
%Dla metody LU\n');

        invA = LU_(A);
        M = A * invA - I;
        delta2 = norm2(M);
        deltaInf = normInf(M);
        delta2norm = norm(M, 2);
        deltaInfnorm = norm(M, inf);
        LUdelta2(k + 1) = delta2;
        LUdeltaInf(k + 1) = deltaInf;
        
%Dla metody inv z MATLABa

        invA = inv(A);
        M = A * invA - I;
        delta2 = norm2(M);
        deltaInf = normInf(M);
        delta2norm = norm(M, 2);
        deltaInfnorm = norm(M, inf);
        invDelta2(k + 1) = delta2;
        invdeltaInf(k + 1) = deltaInf;

      end
    
%wykresy

      figure;
      loglog(zakresX, LLdelta2, zakresX, LUdelta2, zakresX, invDelta2);       
      xlim([zakresX(1), zakresX(end)]);
      legend('\delta2(LL)', '\delta2(LU)', '\delta2(inv)', 'Location', 'NorthWest');
      grid on;
      xlabel('x');
      ylabel('wartoœæ wskaŸnika b³êdu');
      title('Porównanie wyznaczenia normy 2 metodami LL, LU i funkcj¹ inv');

      figure;
      loglog(zakresX, LLdeltaInf, zakresX, LUdeltaInf, zakresX, invdeltaInf);
      xlim([zakresX(1), zakresX(end)]);
      legend('\delta{\infty}(LL)', '\delta{\infty}(LU)', '\delta{\infty}(inv)', 'Location', 'NorthWest');
      grid on;
      xlabel('x');
      ylabel('wartoœæ wskaŸnika b³êdu');
      title('Porównanie wyznaczenia normy {\infty} metodami LL, LU i funkcj¹ inv');
    end
end