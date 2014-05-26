function makeIntralipidPlots
  clear; close all; rng(1);
  addpath(genpath(pwd))

  muAlpha = 0;    % Note, if we know true value than problem is better
  
  numIntralipids = 6;
  numPix = 512;
  
  muFit = zeros(numPix, numIntralipids);
  dataCases = 11:16;
  
  for i = 1:numel(dataCases)
    [I, z, dx, z0, zR, muAlpha, muBeta, muL0, trueMu ] = loadOctData( dataCases(i), false );
    
    I = I(:,150:201);
    I = mean( I, 2 );
    
    mask = findNonZeroMus(I);
    noiseLevel = median( I(mask==0) );
    I = I - noiseLevel;
    I = max( I, 0 );
    I = I.*mask;
    I = I ./ 1000;
    
     eta = 1d-3;
     muFit(:,i) = muFitCVX( I, mask, z, z0, zR, eta );
  end
  
  plot(muFit),ylim([0, 4.5]),
  legend('1.25%', '2.5%', '5%', '10%', '15%', '20%')
  
end