function X = MCP_prox_tnn(Y,rho,gama,c)

[n1,n2,n3] = size(Y);
X = zeros(n1,n2,n3);
Y = fft(Y,[],3);
        
% first frontal slice
[U,S,V] = svd(Y(:,:,1),'econ');
S = diag(S);
X(:,:,1) = U*diag(prox_MCP(S,rho,gama,c))*V';

% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2 : halfn3
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    X(:,:,i) = U*diag(prox_MCP(S,rho,gama,c))*V';
    X(:,:,n3+2-i) = conj(X(:,:,i));
end

% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    X(:,:,i) = U*diag(prox_MCP(S,rho,gama,c))*V';
end
X = ifft(X,[],3);