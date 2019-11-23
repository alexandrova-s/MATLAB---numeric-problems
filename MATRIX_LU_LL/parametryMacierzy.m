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