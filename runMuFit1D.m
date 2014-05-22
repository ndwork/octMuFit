
function runMuFit1D
  clear; close all; rng(1);
  addpath(genpath(pwd))

  muAlpha = 0;    % Note, if we know true value than problem is better

  dataCase = 3;
  %[I, z, dx, z0, zR, muAlpha, muBeta, muL0, trueMu ] = loadOctData( dataCase, false );
load 'data.mat';

  if ~isvector(I)
    if dataCase==0
      trueMu = trueMu(:,1);
      I = I(:,1);
    else
      I = I(:,156:256);
      I = mean( I, 2 );
    end
  end

  if dataCase == 0
    mask = ones( numel(I), 1 );
  else
    mask = findNonZeroMus(I);
  end

  noiseLevel = median( I(mask==0) );
  I = I - noiseLevel;
  I = max( I, 0 );

  eta = 1d5;
  I = I.*mask;
  %muStar = muFitCVX( I, mask, z, z0, zR, eta );
load 'muStar.mat';
  %muFit = muStar;
  tic
  %[muFit, fos, relFos] = muFit1D_CP( I, mask, z, z0, zR, eta, muStar );
  %[muFit, fos, relFos] = muFit1D_CP2( I, mask, z, z0, zR, eta, muStar );
  %[muFit, fos, relFos] = muFit1D_LADMM(I, mask, z, z0, zR, eta, muStar );
  [muFit, fos, relFos] = muFit1D_ADMM(I, mask, z, z0, zR, eta, muStar );
  timeTaken = toc;
  disp(['Time taken: ', num2str(timeTaken)]);
  if exist('fos')
    figure, semilogy(fos);  title('fos');
  end
  if exist('relFos')
    figure, semilogy(relFos, 'LineWidth', 2);
    xlabel('Iteration', 'FontSize', 14)
    ylabel('Relative Error', 'FontWeight', 'bold', 'FontSize', 14)
  end



  figure;
  plot( z, muFit, 'r', 'LineWidth', 2 );
  hold on;  
  if numel( trueMu ) > 0
    plot( z, trueMu, 'b', 'LineWidth', 2 );
    legend( 'CVX Fit', 'True \mu', 'Location', 'NorthWest' );
  end
  ylim([0,4.500]);
  xlim([z(1),z(end)]);
  xlabel('Depth (mm)','FontSize',14);
  ylabel('Attenuation (mm^{-1})', 'FontSize', 14, 'FontWeight', 'bold');

  gamFit = 1 ./ muFit;
  figure;
  plot( z, gamFit, 'r', 'LineWidth', 2 );
  hold on;
  if numel( trueMu ) > 0
    trueGam = 1 ./ trueMu;
    plot( z, trueGam, 'b', 'LineWidth', 2 );
    legend( 'CVX Fit', 'True \gamma' );
  end
  ylim([0 2.5]);
  xlabel('Depth (mm)');
  ylabel('Gamma (mm)');

  figure;
  plot( z, I );
  fitI = mu2I( muFit, z, z0, zR, muAlpha, muBeta, muL0 );
  scale = mean( I(1:50) ./ fitI(1:50) );
  fitI = fitI * scale;
  hold on;
  plot( z, fitI, 'r', 'LineWidth', 2 );
  if numel( trueMu )
    trueI = mu2I( trueMu, z, z0, zR, muAlpha, muBeta, muL0 );
    plot( z, trueI, 'k', 'LineWidth', 2 );
    legend('data','fit','actual');
  end
  xlabel('Depth (mm)');
  ylabel('Intensity');
  
  zless1mm = z < z(1)+0.75; %Find depths for which we are less than 0.75mm away from top surface 
  diffI = norm(I(zless1mm) - fitI(zless1mm));
  disp(['Norm difference between I and fit of I: ' num2str(diffI)])


  dz = z(2)-z(1);
  gam = 1 ./ muFit;
  m = numel(I);
  ut = triu( ones(m-1,m-1) );
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  dz = z(2:end) - z(1:end-1);
  dz = [ dz; dz(end) ];
  fitf = 1/2 * I(1:m-1) .* gam(1:m-1) ./ dz(1:m-1) ...
         - 1/2 * I(m) * gam(m) / dz(m) ...
         + 1/2 * ut * ( ...
           I(1:m-1) ./ dz(1:m-1) .* gam(1:m-1) .* log( g(2:m)./g(1:m-1) ) ...
         );
  f = ut * I(1:m-1);
  figure;
  plot( z(1:m-1), f );
  hold on;
  plot( z(1:m-1), fitf, 'r', 'LineWidth', 2 );
  if numel( trueMu )
    trueI = mu2I( trueMu, z, z0, zR, muAlpha, muBeta, muL0 );
    truef = flipud(cumsum(flipud(trueI)));
    plot( z, truef, 'k', 'LineWidth', 2 );
    legend('data (b)','fit (A\gamma)','actual (A*true\gamma)');
  else
    legend('data (b)','fit (A\gamma)');
  end
  xlabel('Depth (mm)');
  ylabel('f');




end

