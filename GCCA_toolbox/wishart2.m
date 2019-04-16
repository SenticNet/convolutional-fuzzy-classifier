function [sig] = wishart2(sig,df,p);

% taken from Johnson, M. (1987). "multivariate stat. simulation", pp 203-204.

T = zeros(p);

for i=1:p
  for j=1:i-1
    T(i,j) = randn;
  end
  T(i,i) = sqrt(g05fff((df-i+1)/2,2,1));
end

U = chol(sig);
L = U';

sig = L*T*T'*U;

