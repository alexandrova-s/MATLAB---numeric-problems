function D = wykresDet(N,a,rozmiar)
 
    for i = 1:rozmiar  
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