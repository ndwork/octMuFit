
function runMuFit2D
  clear; close all; rng(1);
  addpath(genpath(pwd))

  muAlpha = 0;    % Note, if we know true value than problem is better

  dataCase = 3;
  [I, z, dx, z0, zR, muAlpha, muBeta, muL0, trueMu ] = loadOctData( dataCase, false );
  
  mask = findNonZeroMus(I);

  noiseLevel = median( I( mask==0 ) );
  I = max( I - noiseLevel, 0 );

  muFit = muFit2D_ADMM(I, mask, z, z0, zR, eta, muStar );

end



