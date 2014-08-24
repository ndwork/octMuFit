function muFit = muFit2D_ver( I, z )

  [M, N] = size( I );
  muFit = zeros( M, N );

  for j=1:N
    line = I(:,j);
    muFit(:,j) = muFitVermeer( line, z );
  end
  muFit = max( muFit, 0 );
  muFit( ~isfinite(muFit) ) = 0;

end

