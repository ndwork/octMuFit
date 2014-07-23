


function muFit = muFit2D_mVer_gBlur( I, z, z0, zR )

  muHat = muFit2D_mVer( I, z, z0, zR );

  hSize = 5;
  sigma = 2.0;
  h = fspecial('gaussian', hSize, sigma);
  %h = h( round( hSize / 2 ), : );
  %muFit = filter( h, 1, muHat, [], 2 );
  muFit = imfilter( muHat, h );

  muFit(1:hSize,:) = muHat(1:hSize,:);
  
end

