
function findZo
  clear; close all; rng(1);
  addpath(genpath(pwd))

  dataCase = 9;
  [I, z, dx, trueZ0, zR, muAlpha, muBeta, muL0, ALA, trueMu ] = ...
    loadOctData( dataCase, false );

  if ~isvector(I)
    if dataCase==0
      trueMu = trueMu(:,1);
      I = I(:,1);
    else
      I = I(:,250:274);
      I = mean( I, 2 );
    end
  end

  if dataCase == 0
    mask = trueMu > 0;
  else
    [mask, noiseLevel] = makeTheMask( I, ALA );
    I = I - noiseLevel;
    I = max( I, 0 );
    I(mask==0) = 0;
    %I = I ./ 1000;
    trustThresh = 500;
  end
  
  minZ0 = min(z)-0.2;
  maxZ0 = max(z)+0.2;
  dz0 = 0.02;

  z0s = minZ0:dz0:maxZ0;
  vars = zeros( size( z0s ) );
  for i = 1:numel(z0s)
    z0 = z0s(i);
    
    h = makeConfocalFunction( z, z0, zR );

    mu = muFitModVermeer( I, z, h );
    
    vars(i) = var( mu(I>trustThresh) );
    
  end
  
  [~, z0Indx] = min( vars );
  estimateZ0 = z0s(z0Indx);
  disp(['True / Estimated z0 is: ', num2str(trueZ0), ...
    ' / ', num2str(estimateZ0)]);

  figure;
  hTrue = makeConfocalFunction( z, trueZ0, zR );
  muTrue = muFitModVermeer( I, z, hTrue );
  plot( z, muTrue, 'k' );
  hold on;
  hEst = makeConfocalFunction( z, estimateZ0, zR );
  muEst = muFitModVermeer( I, z, hEst );
  plot( z, muEst, 'r' );
  ylim([0 8]);
  
end
