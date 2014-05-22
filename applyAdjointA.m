
function out = applyAdjointA( y, I, z, dz, z0, zR )

  s = size(I);
  ndimsI = sum( s > 1 );

  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );

  switch ndimsI
    case 1 
      out = applyAdjA1D( y, I, dz, g );
    case 2
      out = applyAdjA2D( y, I, dz, g );
    case 3
      out = applyAdjA3D( y, I, dz, g );
    otherwise
      error('improper use of applyA');
  end

end



function out = applyAdjA3D( y, I, dz, g )
  [M N K] = size( I );
  out = zeros([M N K]);
  for i=1:N
    for j=1:K
      out(:,i,j) = applyAdjA1D( y, I(:,i,j), dz, g );
    end
  end
end

function out = applyAdjA2D( y, I, dz, g )
  [M N] = size(I);
  out = zeros( M, N );
  for i=1:N
    out(:,i) = applyAdjA1D( y(:,i), I(:,i), dz, g );
  end
end

function out = applyAdjA1D( y, I, dz, g )
  M = numel(I);
  out = zeros(M,1);

  Lg = log( g(2:end) ./ g(1:end-1) );
  cumsumY = cumsum(y);
  
  out(1:M-1) = I(1:M-1)/(2*dz) .* Lg .* cumsumY(1:M-1) + ...
    I(1:M-1)/(2*dz) .* y(1:M-1);
  
  out(M) = -I(M)/(2*dz)*cumsumY(M-1);
  
end





