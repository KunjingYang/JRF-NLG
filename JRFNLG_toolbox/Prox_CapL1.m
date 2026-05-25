function U = Prox_CapL1(Z,lambda,v)
U = Z;
w = vecnorm(U);
c = w <= v+lambda/(2*v);
Uc = U(:,c);
wc = w(c);
U(:,c) = Uc.*max(wc-lambda/v,0)./wc;