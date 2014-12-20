
function muFit = muFit2D_whTV( I, z, z0, zR, mask )

  %eta = 1.0;
  %eta = 0.3;  % Best one so far (before divide by Nx)
  eta = 300;  % Best one so far (after dividing by Nx)
  %epsilon = 1d-2;
  epsilon = 1d-1;

  [M, N] = size( I );

  muHat = muFit2D_mVer( I, z, z0, zR );

  weight = 1 ./ ( I + epsilon );

  Dx = -eye(N);
  for i=1:N-1
    Dx(i,i+1) = 1;
  end
  Dx(N,:) = 0;

  muFit = zeros( M, N );
  for i=1:M
    if mod(i,50)==0
      disp(['Working on row ', num2str(i), ' of ', num2str(M)]);
    end
    
    muHatRow = muHat(i,:)';
    weightRow = weight(i,:)';
    maskRow = mask(i,:)';

    if max( maskRow ) == 0
      muFit(i,:) = 0;
    else
    
      Nx = sum(maskRow);

      cvx_begin quiet
        variable muRow(N,1)
        minimize 0.5 * sum(( muRow - muHatRow ).^2) + ...
          eta/Nx * norm( weightRow .* maskRow .* ( Dx * muRow ), 1 )
        %minimize 0.5 * norm( muRow - muHatRow, 1 ) + ...
        %  eta/Nx * norm( weightRow .* maskRow .* ( Dx * muRow ), 1 )
        subject to
          muRow >= 0
      cvx_end

      muFit(i,:) = muRow;
      
    end
  end

end

