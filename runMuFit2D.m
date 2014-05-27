
function runMuFit2D
  clear; close all; rng(1);
  addpath(genpath(pwd));

  muAlpha = 0;    % Note, if we know true value than problem is better

  dataCase = 0;
  [I, z, dx, z0, zR, muAlpha, muBeta, muL0, trueMu ] = loadOctData( dataCase, false );
%load 'data.mat';

  %I = I(:,150:200);

  mask = findNonZeroMus(I);

  noiseLevel = median( I( mask==0 ) );
  I = max( I - noiseLevel, 0 );
  I = I .* mask;
  I = I ./ 1000;

  etaz = 1d-3;
  etax = 1d-3;
  profile clear;
  profile on;
  tic;
  [muFit, diagnostics ]= muFit2D_ADMM(I, mask, z, dx, z0, zR, etaz, etax );
  timeTaken = toc;
  profile off;
  disp(['Time taken (s):', num2str(timeTaken)]);

  figure, imshow( muFit, [0 4.5] );
  figure, semilogy( fos ); title('fos');  xlabel('ADMM Iteration');
  figure, plot( muFit(:,2) );  title('muFit col 25'); ylim([0 4.5]);

  figure, plot( cgNIters );  title('CG N Iterations');  xlabel('ADMM Iteration');
  figure, plot( cgRelErrors );  title('CG Relative Errors');  xlabel('ADMM Iteration');

  profile viewer;
end



