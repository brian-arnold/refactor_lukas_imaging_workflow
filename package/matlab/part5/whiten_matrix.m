function M = whiten_matrix(M, factor)

M = bsxfun(@minus, M, mean(M));
A = M'*M;

[V,D] = eig(A);
M = M*V*diag(1./(diag(D)+factor).^(1/2))*V';

