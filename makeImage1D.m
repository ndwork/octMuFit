function makeImage1D
  clear; close all; rng(1);
  addpath(genpath(pwd))


  dataCase = 10;
  [I, z, ~, z0, zR, ~, ~, ~, trueMu ] = loadOctData( dataCase, false );

  if dataCase == 7
    I = I(100:end,:);
  end

  if dataCase == 0
    mask = ones( numel(I), 1 );
    eta = 1d4
  else
    mask = findNonZeroMus(I);
    noiseLevel = median( I(mask==0) );
    I = I - noiseLevel;
    I = max( I, 0 );
    I = I.*mask;
    I = I ./ 1000;
    eta = 1;
  end

  muFit = zeros(size(I));
  [m, n] = size(I);


  numColAvg = 15;
  halfNCol = floor(numColAvg/2);

  % Horizontally average the data
  h = fspecial('gaussian', [1 numColAvg], 2);
  %h = fspecial('average', [1 numColAvg]);
  Iavg = imfilter( I, h );

  for i = 1:n
    if(mod(i, 1) == 0)
      disp(['Working on column ', num2str(i)]);
    end

    % Average column of interest with adjacent columns
    if( i <= halfNCol || i >= (n-halfNCol) ) continue; end;

    %Iavg = mean(I(:,i-halfNCol:i+halfNCol), 2);
    muFit(:,i) = muFitCVX( Iavg(:,i), mask(:,i), z, z0, zR, eta );
  end

  imshow(muFit, [0, 4.5])
end

