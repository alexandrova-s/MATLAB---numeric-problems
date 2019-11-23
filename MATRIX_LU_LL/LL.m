function LL = LL(A)
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