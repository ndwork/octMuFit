
function muFit = muFit2D_DRC( I, z, z0, zR )

  [M, N] = size( I );
  muFit = zeros( M, N );

  h = makeConfocalFunction( z, z0, zR );

  for j=1:N
    if mod(j,5)==0
      disp(['Working on line ', num2str(j), ' of ', num2str(N)]);
    end
    line = I(:,j);
    muFit(:,j) = muFitDRC( line, z, h );
  end

  muFit = max( muFit, 0 );
end
