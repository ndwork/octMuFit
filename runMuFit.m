
function runMuFit
  clear; close all; rng(1);
  addpath(genpath(pwd))

  alpha = 0;    % Note, if we know true value than problem is better

  dataCase = 0;
  [I, z, dx, z0, zR, alpha, beta, L0, trueMu ] = loadOctData( dataCase, false );

  %mask = findNonZeroMus( I );
  mask = ones( size(I) );

  dz = z(2) - z(1);
  muFit = muFit2D_CP( I, mask, z, dz, dx, z0, zR );
  %muFit = muFit2D_ADMM( I, mask, z, dz, dx, z0, zR );

  colOfInterest = 50;
  I = I(:,colOfInterest);
  muFit = muFit(:,colOfInterest);
  trueMu = trueMu(:,colOfInterest);

  z_mm = z * 1000;
  muFit_mm = muFit / 1000;
  trueMu_mm = trueMu / 1000;


  figure;
  plot( z_mm, muFit_mm, 'r', 'LineWidth', 2 );
  hold on;  
  if numel( trueMu ) > 0
    plot( z_mm, trueMu_mm, 'b', 'LineWidth', 2 );
    legend( 'CP Fit', 'True \mu', 'Location', 'NorthWest' );
  end
  ylim([0,4.500]);
  xlim([z_mm(1),z_mm(end)]);
  xlabel('Depth (mm)');
  ylabel('Attenuation (mm^{-1})');

  figure;
  plot( z_mm, I );
  fitI = mu2I( muFit, z, z0, zR, alpha, beta, L0 );
  scale = mean( I(1:50) ./ fitI(1:50) );
  fitI = fitI * scale;
  hold on;
  plot( z_mm, fitI, 'r', 'LineWidth', 2 );
  if numel( trueMu )
    trueI = mu2I( trueMu, z, z0, zR, alpha, beta, L0 );
    plot( z_mm, trueI, 'k', 'LineWidth', 2 );
    legend('data','fit','actual');
  end
  xlabel('Depth (mm)');
  ylabel('I');
  
end
