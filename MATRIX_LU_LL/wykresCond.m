function C = wykresCond(N,a,size)
 
    for i =1:size  
        P = parametryMacierzy(N(i),a);
        q=P(1,:);
        e=P(3,:);
        semilogy(q,e);
        hold on ;
    end
    grid on;
    
    xlabel('Parametr alfa');
    ylabel('Wskaünik uwarunkowania macierzy');
    title('Wykres cond od alfa');
    
    hold off;
end