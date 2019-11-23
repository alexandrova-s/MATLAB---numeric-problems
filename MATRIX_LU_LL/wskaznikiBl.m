function Z5 = wskaznikiBl(N, rozmiar)

    K = 22;
    
    zakresK = 0 : (K - 1);
    zakresX = zakresK / 300
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

        invA = LL(A);
        M = A * invA - I;
        delta2 = norm2(M);
        deltaInf = normInf(M);
        delta2norm = norm(M, 2);
        deltaInfnorm = norm(M, inf);
        LLdelta2(k + 1) = delta2;
        LLdeltaInf(k + 1) = deltaInf;
        
%Dla metody LU\n');

        invA = LU(A);
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
      semilogy(zakresX, LLdelta2, zakresX, LUdelta2, zakresX, invDelta2);       
      xlim([zakresX(1), zakresX(end)]);
      legend('\delta2(LL)', '\delta2(LU)', '\delta2(inv)', 'Location', 'NorthWest');
      grid on;
      xlabel('x');
      ylabel('wartoœæ wskaŸnika b³êdu');
      title('Porównanie wyznaczenia normy 2 metodami LL, LU i funkcj¹ inv');

      figure;
      semilogy(zakresX, LLdeltaInf, zakresX, LUdeltaInf, zakresX, invdeltaInf);
      xlim([zakresX(1), zakresX(end)]);
      legend('\delta{\infty}(LL)', '\delta{\infty}(LU)', '\delta{\infty}(inv)', 'Location', 'NorthWest');
      grid on;
      xlabel('x');
      ylabel('wartoœæ wskaŸnika b³êdu');
      title('Porównanie wyznaczenia normy {\infty} metodami LL, LU i funkcj¹ inv');
end