
function runMuFit2D
  clear; close all; rng(1);
  addpath(genpath(pwd));

  muAlpha = 0;    % Note, if we know true value than problem is better

  dataCase = 5;
  [I, z, dx, z0, zR, muAlpha, muBeta, muL0, ALA, trueMu ] = ...
    loadOctData( dataCase, false );

  %I = I(:,150:200);

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
  %muFit = muFit2D_TV( I, z, z0, zR );
  %muFit = muFit2D_whTV( I, z, z0, zR, mask );
  muFit = weightedTvDenoise_CP( muFit_mVer, paramsCP );
  %muFit = muFit2D_vReg( I, z, z0, zR, mask );
  %muFit = muFit2D_mVer_gBlur( I, z, z0, zR );
  timeTaken = toc;
  %profile off;
  disp(['Time taken (s):', num2str(timeTaken)]);

  figure, imshow( muFit, [0 5.0] );

  if exists( 'diagnostics' )
    figure, semilogy( fos ); title('fos');  xlabel('ADMM Iteration');
    figure, plot( muFit(:,2) );  title('muFit col 25'); ylim([0 4.5]);

    figure, plot( cgNIters );  title('CG N Iterations');  xlabel('ADMM Iteration');
    figure, plot( cgRelErrors );  title('CG Relative Errors');  xlabel('ADMM Iteration');
  end

  %profile viewer;
end

