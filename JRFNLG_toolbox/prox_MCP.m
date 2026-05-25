function y = prox_MCP(A,mu,gama,c)
   % 注意这里gama要大于rho
   y = sign(A).*(abs(A)-c*mu)/(1-mu/gama);
   C = abs(A)-mu;
   C1 = C<=0;
   y(C1) = 0;
   D = abs(A)-c*gama;
   D1 = D>=0;
   y(D1) = A(D1);
