function LU = LU(A)
n = size(A);
[L, U, P] = lu(A);
I=eye(n);
Y=ones(n);
X=ones(n);
for k = 1:n
    for i = 1:n 
        sumL=0;        
        for j=1:i-1        
            if(i>1)
            sumL=sumL+L(i,j)*X(j,k);
            end             
        end        
        X(i,k)=(I(i,k)-sumL)/L(i,i);        
    end    
end

for a = 1:n
    for b = n:-1:1     
        sumU=0;        
        for c = n:-1:b+1
           sumU=sumU+U(b,c)*Y(c,a);
        end        
        Y(b,a)=(X(b,a)- sumU)/U(b,b);       
    end    
end
LU = Y*P;

end