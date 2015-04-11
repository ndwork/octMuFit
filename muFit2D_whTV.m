
function muFit = muFit2D_whTV( I, z, z0, zR, mask )

  eta = 3;    % good value
  epsilon = 1d-3;

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

      wmRow = muHatRow .* weightRow;
      cvx_begin quiet
        variable muRow(N,1)
        minimize 0.5 * sum(( muRow - muHatRow ).^2) + ...
          eta/Nx * norm( wmRow .* ( Dx * muRow ), 1 )
        %minimize 0.5 * norm( muRow - muHatRow, 1 ) + ...
        %  eta/Nx * norm( wmRow .* ( Dx * muRow ), 1 )
        subject to
          muRow >= 0
      cvx_end

      muFit(i,:) = muRow;
      
    end
  end

end

