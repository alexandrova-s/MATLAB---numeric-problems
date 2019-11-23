function delta = norm2(A)

    delta = max(sqrt(eigs(A' * A)));

end