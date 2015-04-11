
function muFit = muFit2D_mVer( I, z, z0, zR )

  [M, N] = size( I );
  muFit = zeros( M, N );

  h = makeConfocalFunction( z, z0, zR );

  for j=1:N
%for j=495:495
%for j=100:100
    line = I(:,j);
    muFit(:,j) = muFitModVermeer( line, z, h );
  end

end

