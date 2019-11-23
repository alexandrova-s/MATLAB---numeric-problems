clear
close all
clc
 
N = [10, 20, 30];
sizeN = 3;
K = 8;

%wzór funkcji - funkcja anonimowa

f = @(x)((x+1/3).^2 + exp(-x-2));

%Zadanie 1 - sporz¹dzanie wykresu funkcji

for i = 1:sizeN
    x = linspace(-1,1,N(i));
    y = f(x);
    
    figure("Name", "Wykres funkcji", "NumberTitle", "off");
    plot(x,y);
    
    hold on;
    grid on;
    
    xlabel('x');
    ylabel('y');
    title('f(x)');    
    
hold off;

end


%Zadanie 2 - aproksymacja funkcji

for i = 1:sizeN
    
    orgX = linspace(-1,1,N(i));
    orgY = f(orgX);
    
    W = macierzFI(N(i), K);
    
    WWt = W'*W;
    Wy = W' * orgY';
    c = WWt\Wy;
    
    aproxX = linspace(-1,1,1000);
    F = macierzFI(1000, K);
    aproxY = F * c;    
    
    figure("Name", "Aproksymowana funkcja", "NumberTitle", "off");
    plot(orgX, orgY, '-*', aproxX, aproxY, '-');
    
    xlabel('x');
    ylabel('y');
    title('Aproksymacja funkcji f(x)'); 
    grid on;
    
end


%Zadanie 3 - wysresy b³êdów

N = 50;

x = linspace(-1,1,1000)';
y = f(x);

Error2 = zeros(N, N-1);
ErrorInf = zeros(N, N-1);

for n = 5:N
    for k=1:n-1
        orgX = linspace(-1,1,n);
        orgY = f(orgX);
        
        WforErr = macierzFI(n,k);
        WWt = WforErr'*WforErr;
        Wy = WforErr' * orgY';
        c = WWt\Wy;
        
        aproxX = linspace(-1,1,1000);
        F = macierzFI(1000, k);
        aproxY = F * c;    
                
        Error2(n,k) = norm(aproxY-y, 2)/norm(y, 2); 
        ErrorInf(n,k) = norm(aproxY-y, 'inf')/norm(y, 'inf'); 
        
    end
end
 
 %wykres b³êdu œredniokwadratowego
 figure("Name", "WskaŸnik b³êdu œredniokwadratowego", "NumberTitle", "off");
 mesh(log10(Error2)) ;   
 xlabel("K"); 
 ylabel("N"); 
 zlabel("wskaŸnik b³êdu œredniokwadratowego"); 
 grid on;     
 
 %wykres b³êdu maksymalnego
 figure("Name", "WskaŸnik b³êdu maksymalnego", "NumberTitle", "off");
 mesh(log10(ErrorInf));
 xlabel("K"); 
 ylabel("N");   
 zlabel("wskaŸnik b³êdu maksymalnego");
 grid on;      

 %Zadanie 4 - aproksymacja przy zaburzonych danych
 
 minNorma2 = [];

for dy = logspace(-5, -1)
    for N = 5:50
        xc = linspace(-1, 1, 1000);
        yc = f(xc);
        yc = yc + randn(size(yc))*dy;

        Error2 = zeros(N, N-1);
        
        for K = 1:N-1        
            xa = linspace(-1, 1,N);
            ya = f(xa);
            
            W = macierzFI(N,K);
            WWt = W'*W;
            Wy = W' * ya';
            c = WWt\Wy;

            F = macierzFI(1000, K);
            aproxY = F * c;   

            Error2(N,K) = norm(aproxY-yc', 2)/norm(yc', 2); 
        end
    end
    
    minNorma2 = [minNorma2;min(min(Error2(Error2>0)))];
end

figure()
plot(linspace(-5,-1,50),minNorma2, '*');
hold on;
p = polyfit(linspace(-5,-1,50),minNorma2',7);
plot(linspace(-5,-1,50),polyval(p,linspace(-5,-1,50)),'r:');
xlabel("\delta{y}"); 
ylabel("\delta2{MIN}"); 
grid on;            



%funkcja generuj¹ca macierz fi

function M = macierzFI(N,K)

    xn = linspace(-1,1,N);
    xk = linspace(-1,1,K);
    M=zeros(N,K);
    
    for k = 1:K
        for n = 1:N           
            s = abs(xn(n)-xk(k));
            if s < 1/2                
                M(n,k) = 16*abs(xn(n)-xk(k))^3 - 12*abs(xn(n)-xk(k))^2 + 1;
            end
            if s >= 1/2
                M(n,k) = 0;
            end 
        end
    end

end
 