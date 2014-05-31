
function wrongZos
  clear; close all; rng(1);
  addpath(genpath(pwd))

  muAlpha = 0;    % Note, if we know true value than problem is better

  dataCase = 7;
  [I, z, dx, z0, zR, muAlpha, muBeta, muL0, trueMu ] = loadOctData( dataCase, false );

  figure;
  hold all;
  linestyleOrder = {'--','-',':'};

  numColAvg = 15;
  h = fspecial('gaussian', [1 numColAvg], 3);
  Iavg = imfilter( I, h );
  mask = findNonZeroMus(Iavg);
  I = Iavg(:,260);
  mask = mask(:,260);


  noiseLevel = median( I(mask==0) );
  I = I - noiseLevel;
  I = max( I, 0 );
  I(mask==0) = 0;
  I = I ./ 1000;
  eta = 1;

  iter = 1;
  dStep = 0.3;
  %for diff = -dStep:dStep:dStep
  for diff = 0:0
    thisz0 = z0 + diff;

    muFit = muFitCVX( I, mask, z, thisz0, zR, eta );
    plot( z, muFit, 'LineStyle', linestyleOrder{iter}, 'LineWidth', 1.5 );
    
    iter = iter+1;
  end

  legend( '\Delta z = -0.4', '\Delta z = 0', '\Delta z = 0.4');

  xlabel('Depth (mm)', 'FontSize', 14);
  ylabel('Attenuation Coefficient (mm^{-1})', 'FontSize', 14);
  axis([0 max(z) 0 8]);


end

