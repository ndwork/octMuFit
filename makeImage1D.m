function makeImage1D
  clear; close all; rng(1);
  addpath(genpath(pwd))


  dataCase = 10;
  [I, z, ~, z0, zR, ~, ~, ~, trueMu ] = loadOctData( dataCase, false );


  numColAvg = 15;
  
  % Horizontally average the data
  h = fspecial('gaussian', [1 numColAvg], 3);
  %h = fspecial('average', [1 numColAvg]);
  Iavg = imfilter( I, h );
  
  if dataCase == 0
    mask = trueMu > 0;
    dMu = trueMu(2:end) - trueMu(1:end-1);
    eta = 1d4;
  else
    mask = findNonZeroMus(Iavg);
    noiseLevel = median( Iavg(mask==0) );
    Iavg = Iavg - noiseLevel;
    Iavg = max( Iavg, 0 );
    Iavg(mask==0) = 0;
    Iavg = Iavg ./ 1000;
    eta = 1d-1;
  end

  muFit = zeros(size(Iavg));
  [m, n] = size(Iavg);

  halfNCol = floor(numColAvg/2);
  for i = halfNCol+1:n-halfNCol-1
    if(mod(i, 1) == 0)
      disp(['Working on column ', num2str(i)]);
    end

    muFit(:,i) = muFitCVX( Iavg(:,i), mask(:,i), z, z0, zR, eta );
  end

  imshow(muFit, [0, 4.5])
end

