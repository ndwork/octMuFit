
function runMuFit2D
  clear; close all; rng(1);
  addpath(genpath(pwd));

  muAlpha = 0;    % Note, if we know true value than problem is better

  dataCase = 13;
  %dataCase = 0;   % Simulation
  [I, z, dx, z0, zR, muAlpha, muBeta, muL0, ALA, trueMu ] = ...
    loadOctData( dataCase, false );

  if dataCase == 0
    mask = ones( size(I) );
    noiseLevel = 0;
  else
    mask = findNonZeroMus(I);
    noiseLevel = median( I( mask==0 ) );
  end

  I = max( I - noiseLevel, 0 );
  I = I .* mask;
  I = I ./ 1000;

  %TV Parameters
  maxIter = 1000;
  tau = .02;  % Actually, the CP codes choose step sizes tau = sigma = 1/nrmK.
  overRelax = 1.9;
  epsilon = 1d-1;
  weight = 1 ./ ( I + epsilon );
  gamma = .2*weight;
  showTrigger = pi;  %by making irrational, never shows
  paramsCP = struct('gamma',gamma,'tau',tau,'maxIter',maxIter,...
            'showTrigger',showTrigger,'theta',1,'overRelax',overRelax);

  etaz = 1d-3;
  etax = 1d-3;
  %profile clear;
  %profile on;
  tic;
  %[muFit, diagnostics ]= muFit2D_ADMM(I, mask, z, dx, z0, zR, etaz, etax );
  muFit_mVer = muFit2D_mVer( I, z, z0, zR );

  %muFit = muFit_mVer;
  %muFit = muFit2D_TV( I, z, z0, zR );
  muFit = muFit2D_whTV( I, z, z0, zR, mask );
  %muFit = weightedTvDenoise_CP( muFit_mVer, paramsCP );
  %muFit = muFit2D_vReg( I, z, z0, zR, mask );
  %muFit = muFit2D_mVer_gBlur( I, z, z0, zR );
  timeTaken = toc;
  %profile off;
  disp(['Time taken (s):', num2str(timeTaken)]);

  if dataCase >= 11 && dataCase <= 16
    dz = z(2)-z(1);
    [~, N] = size(I);
    midCol = floor(N/2);
    skinLoc = find( mask(:,midCol)==1, 1 );
    depth5mm = skinLoc + floor(0.5/dz);
    meanI = mean(I,2);
    muFit_faber = muFitFaber( meanI, [skinLoc,depth5mm], z, z0, zR );
    %plot( muFit_faber, 'r', 'LineWidth', 2 );
    %hold on;
    %plot( muFit_mVer(:,midCol), 'b', 'LineWidth', 1 );  ylim([0 5]);
    meanFaber = mean( muFit_faber(skinLoc:depth5mm) );
    meanFit = mean( muFit(skinLoc:depth5mm,midCol) );
    disp(['MuFit Faber: ', num2str( meanFaber ) ]);
    disp(['MuFit mVer: ', num2str( meanFit ) ]);
  end

  figure, imshow( muFit, [0 5.0] );

  if exist( 'diagnostics' )
    figure, semilogy( fos ); title('fos');  xlabel('ADMM Iteration');
    figure, plot( muFit(:,2) );  title('muFit col 25'); ylim([0 4.5]);

    figure, plot( cgNIters );  title('CG N Iterations');  xlabel('ADMM Iteration');
    figure, plot( cgRelErrors );  title('CG Relative Errors');  xlabel('ADMM Iteration');
  end

  %profile viewer;
end

