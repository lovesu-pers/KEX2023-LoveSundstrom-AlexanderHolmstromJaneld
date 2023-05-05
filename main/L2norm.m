function L2 = L2norm(v1,v2)
    N = length(v1);
    L2 = sqrt(1/N * sum( abs(v1-v2).^2 ) );
end