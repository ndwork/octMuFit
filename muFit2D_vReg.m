
function muFit = muFit2D_vReg( I, z, z0, zR, mask )

  %eta = 1.0;
  %eta = 0.3;  % Best one so far (before divide by Nx)
  eta = 100;
  %epsilon = 1d-2;
  epsilon = 1d-1;

  [M, N] = size( I );
  muHat = zeros( M, N );

  h = makeConfocalFunction( z, z0, zR );
  for j=1:N
    muHat(:,j) = muFitModVermeer( I(:,j), z, f, h );
  end
  muHat = max( muHat, 0 );
  muHat( ~isfinite(muHat) ) = 0;

  weight = 1 ./ ( I + epsilon );

  Dx = -eye(M);
  for i=1:M-1
    Dx(i,i+1) = 1;
  end
  Dx(M,:) = 0;

  muFit = zeros( M, N );
  for j=1:N
    disp(['Working on col ', num2str(j), ' of ', num2str(N)]);
    
    muHatCol = muHat(:,j);
    weightCol = weight(:,j);
    maskCol = mask(:,j);

    Nx = sum( maskCol );
    
    cvx_begin quiet
      variable muCol(M,1)
      minimize 0.5 * sum(( muCol - muHatCol ).^2) + ...
        eta/Nx * norm( weightCol .* maskCol .* ( Dx * muCol ), 1 )
      subject to
        muCol >= 0
    cvx_end

    muFit(:,j) = muCol;
  end

end

