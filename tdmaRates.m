function R = tdmaRates(NC,NU,h,alpha,Pmax,nvar)
R = zeros(NU,NC);
g = abs(h).^2;
alpha = ones(NU, NC)/NU;
for c=1:NC
    for u=1:NU
        SINR = (g(u,c,c)*Pmax)/(sum(g(u,c,:).*Pmax)+nvar);
        R(u,c) = alpha(u,c)*log2(1+SINR);
    end
end