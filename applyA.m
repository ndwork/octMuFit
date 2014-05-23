

function out = applyA( gamma, I, dz, z, z0, zR )
  % I is the intensity data
  % dz is a scalar representing the size of each pixel in the z (depth)
  %   direction
  % gamma the inverse attenuation coefficient matrix, the same size as I

  s = size(I);
  ndimsI = sum( s > 1 );

  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );

  switch ndimsI
    case 1
      M = numel(I);
      ut = triu( ones(M-1,M-1) );
      out = applyA1D( gamma, I, dz, g, ut );
    case 2
      out = applyA2D( gamma, I, dz, g );
    case 3
      out = applyA3D( gamma, I, dz, g );
    otherwise
      error('improper use of applyA');
  end

end


function out = applyA3D( gamma, I, dz, g )
  [M N K] = size( I );
  out = zeros([M N K]);
  
  ut = triu( ones(M-1,M-1) );

  tmp = repmat( log( g(2:M)./g(1:M-1) ), [1 N, K] );
  
  out(1:M-1,:,:) = 1/2 * I(1:M-1,:) .* gamma(1:M-1,:) ./ dz ...
    - 1/2 * I(M) * gamma(M) ./ dz ...
    + 1/2 * ut * ( ...
      I(1:M-1,:) ./ dz ./ gamma(1:M-1,:) .* tmp ...
    );
end


function out = applyA2D( gamma, I, dz, g )
  [M, N] = size(I);
  out = zeros( M, N );
  
  ut = triu( ones(M-1,M-1) );

  tmp = repmat( log( g(2:M)./g(1:M-1) ), [1 N] );
  
  out(1:M-1,:) = 1/2 * I(1:M-1,:) .* gamma(1:M-1,:) ./ dz ...
    - 1/2 * I(M) * gamma(M) ./ dz ...
    + 1/2 * ut * ( ...
      I(1:M-1,:) ./ dz ./ gamma(1:M-1,:) .* tmp ...
    );
end


function out = applyA1D( gamma, I, dz, g, ut )

  M = numel(I);

  out = zeros( M, 1 );
  out(1:M-1) = 1/2 * I(1:M-1) .* gamma(1:M-1) ./ dz ...
    - 1/2 * I(M) * gamma(M) / dz ...
    + 1/2 * ut * ( ...
      I(1:M-1) ./ dz .* gamma(1:M-1) .* log( g(2:M)./g(1:M-1) ) ...
    );

end

