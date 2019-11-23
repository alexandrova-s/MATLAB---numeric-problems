%Generacja Macierzy 
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

