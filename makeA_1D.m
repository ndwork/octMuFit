
function A = makeA_1D( I, dz, g )

  M = numel(I);

  A = zeros(M,M);

  Lg = log( g(2:M) ./ g(1:M-1) );

  for i=1:M-1
    A(i,i) = 1/(2*dz) * I(i) * ( 1 + Lg(i) );
    A(i,M) = -1/(2*dz) * I(M);
  end

  for i=1:M-1
    for j=i+1:M-1
      A(i,j) = 1/(2*dz) * I(j) * Lg(j);
    end
  end

end
