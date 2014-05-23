
function runMuFit2D
  clear; close all; rng(1);
  addpath(genpath(pwd));

  muAlpha = 0;    % Note, if we know true value than problem is better

  dataCase = 3;
  [I, z, dx, z0, zR, muAlpha, muBeta, muL0, trueMu ] = loadOctData( dataCase, false );

  I = I(:,150:200);

  mask = findNonZeroMus(I);

  noiseLevel = median( I( mask==0 ) );
  I = max( I - noiseLevel, 0 );
  I = I .* mask;

  etaz = 1d5;
  etax = 1d7;
  profile clear;
  profile on;
  tic;
  [muFit, fos] = muFit2D_ADMM(I, mask, z, dx, z0, zR, etaz, etax );
  timeTaken = toc;
  profile off;
  disp(['Time taken (s):', num2str(timeTaken)]);

  figure, imshow( muFit, [0 4.5] );
  figure, plot( fos ); title('fos');  
  figure, plot( muFit(:,2) );  title('muFit col 25'); ylim([0 4.5]);
  
  profile viewer;
end



