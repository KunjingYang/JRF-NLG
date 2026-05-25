function y = prox_MCP1(gama,A,miu)
   % 注意这里gama要大于1 
   y = sign(A).*(abs(A)-miu)/(1-1/gama);
   C = abs(A)-miu;
   C1 = C<=0;
   y(C1) = 0;
   D = abs(A)-gama*miu;
   D1 = D>=0;
   y(D1) = A(D1);
