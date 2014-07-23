
function muFit = muFit2D_mVer( I, z, z0, zR )

  [M, N] = size( I );
  muFit = zeros( M, N );

  h = makeConfocalFunction( z, z0, zR );
  for j=1:N
    line = I(:,j);
    muFit(:,j) = muFitModVermeer( line, z, h );
  end
  muFit = max( muFit, 0 );
  muFit( ~isfinite(muFit) ) = 0;
  
end

