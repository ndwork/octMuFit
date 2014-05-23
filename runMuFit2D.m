
function runMuFit2D
  clear; close all; rng(1);
  addpath(genpath(pwd))

  muAlpha = 0;    % Note, if we know true value than problem is better

  dataCase = 3;
  [I, z, dx, z0, zR, muAlpha, muBeta, muL0, trueMu ] = loadOctData( dataCase, false );
  
  I = I(:,150:200);
  
  mask = findNonZeroMus(I);

  I = I .* mask;
  noiseLevel = median( I( mask==0 ) );
  I = max( I - noiseLevel, 0 );

  eta = 1d5;
  [muFit fos] = muFit2D_ADMM(I, mask, z, dx, z0, zR, eta );

  figure, imshow( muFit, [0 4.5] );
  figure, plot( fos ); title('fos');
  
end



