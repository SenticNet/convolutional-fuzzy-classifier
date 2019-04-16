function [sig] = wishart(sig, df, k)

% WISHART FOR k X k matrix, degrees of freedom = df

Ep = zeros(k);

for i=1:k
  for j=i:k
    Ep(i,j) = randn;
  end
end

V = zeros(k);

for i=1:k

  Pittilon = g05fff((df - i + 1)/2, 2, 1);
  V(i,i) = Pittilon;

  for l=1:i

    V(i,i) = V(i,i) + (Ep(l,i)^2);
 
  end

  for j=i+1:k

    V(i,j) = Ep(i,j) * sqrt(Pittilon);

    for m=1:i

       V(i,j) = V(i,j) + (Ep(m,i)*Ep(m,j));

    end

    V(j,i) = V(i,j);

  end

end

U = chol(sig);
L = U';

Ep = L*V;

sig = Ep*U;



