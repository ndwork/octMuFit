
function runMuFit2D
  clear; close all; rng(1);
  addpath(genpath(pwd));

  muAlpha = 0;    % Note, if we know true value than problem is better

  dataCase = 16;
  %dataCase = 0;   % Simulation
  [I, z, dx, z0, zR, muAlpha, muBeta, muL0, lambda, deltaLambda, ...
    dLambda, ALA, trueMu ] = loadOctData( dataCase, false );

  if dataCase == 0
    mask = ones( size(I) );
    noiseLevel = 0;
  else
    mask = findNonZeroMus(I);
    noiseLevel = median( I( mask==0 ) );
  end

  %I = I - noiseLevel;
  %I = max( I, 0 );
  I = I .* mask;

  %profile clear;
  %profile on;
  tic;
  %muFit = muFit2D_ver( I, z );
  muFit = muFit2D_mVer( I, z, z0, zR );
  %muFit = muFit2D_DRC( I, z, z0, zR );
  %muFit = muFit2D_mVer_gBlur( I, z, z0, zR );
  %muFit = muFit2D_TV( I, z, z0, zR );
  %muFit = muFit2D_whTV( I, z, z0, zR, mask );
  %muFit = muFit2D_vReg( I, z, z0, zR, mask );
  timeTaken = toc;
  %profile off;
  disp(['Time taken (s):', num2str(timeTaken)]);

  if dataCase >= 11 && dataCase <= 16
    avgDepth = 0.25;
    dz = z(2)-z(1);
    [~, N] = size(I);
    midCol = floor(N/2);
    skinLocs = zeros(1,N);
    for i=1:N, skinLocs(i) = find( mask(:,i)==1, 1 ); end;
    depthMM = skinLocs + floor(avgDepth/dz);
    meanI = mean(I,2);
    midSkinLoc = skinLocs(midCol);
    muFit_faber = muFitFaber( meanI, [midSkinLoc], z, z0, zR );
    %plot( muFit_faber, 'r', 'LineWidth', 2 );
    %hold on;
    %plot( muFit_mVer(:,midCol), 'b', 'LineWidth', 1 );  ylim([0 5]);
    meanFaber = mean( muFit_faber(midSkinLoc:depthMM) );
    meanFit = mean( muFit(midSkinLoc:depthMM,midCol) );
    disp(['MuFit Faber: ', num2str( meanFaber ) ]);
    disp(['MuFit: ', num2str( meanFit ) ]);
    medMuFit = median( muFit, 2 );
    medFit = mean( medMuFit(midSkinLoc:depthMM) );
    disp(['MuFit Median: ', num2str( medFit ) ]);

    avgMus = muFit(skinLocs:depthMM,:);
    avgMus = mean( avgMus, 1 );
    disp(['Avg Mu for all cols: ', num2str(mean(avgMus)) ]);
    disp(['Med Mu for all cols: ', num2str(median(avgMus)) ]);
    disp(['StdDev for all cols: ', num2str(std(avgMus)) ]);
    disp(['Max Mu for all cols: ', num2str(max(avgMus)) ]);
    disp(['Min Mu for all cols: ', num2str(min(avgMus)) ]);
  end

  figure, imshow( muFit, [0 3.0] );

  if dataCase == 0
    figure; plot( z, trueMu(:,50), 'k' ); ylim([0 5]);
    hold on; plot( z, muFit(:,50), 'b' );
    [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit,trueMu,z);
    disp(['Mean ETB Depth: ', num2str(meanEtbDepth)]);
    disp(['Mean V Depth: ', num2str(meanVDepth)]);
  end

  if exist( 'diagnostics', 'var' )
    figure, semilogy( fos ); title('fos');  xlabel('ADMM Iteration');
    figure, plot( muFit(:,2) );  title('muFit col 25'); ylim([0 4.5]);

    figure, plot( cgNIters );  title('CG N Iterations');  xlabel('ADMM Iteration');
    figure, plot( cgRelErrors );  title('CG Relative Errors');  xlabel('ADMM Iteration');
  end

  %profile viewer;
end

