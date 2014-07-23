
function simMuFit2D( )
  clear; close all; rng(100);
  addpath(genpath(pwd));

  dataCase = 0;  % dataCase=0 for simulation
  [I, z, dx, z0, zR, muAlpha, muBeta, muL0, ALA, trueMu ] = ...
    loadOctData( dataCase, false );

  mask = ones( size(I) );
  noiseLevel = 0;

  [M N] = size( I );
  %I = max( I - noiseLevel, 0 );
  %I = I .* mask;
  %I = I ./ 1000;  % don't divide by 1000 for simulation

  col = round(N/2);
  lineI = I(:,col);
  mask = trueMu > 0;
  lineMask = mask(:,col);

  %TV Parameters
  maxIter = 1000;
  tau = .02;  % Actually, the CP codes choose step sizes tau = sigma = 1/nrmK.
  overRelax = 1.9;
  epsilon = 1d-1;
  weight = 1 ./ ( I + epsilon );
  gamma = .5*weight;
  showTrigger = pi;  %by making irrational, never shows
  paramsCP = struct('gamma',gamma,'tau',tau,'maxIter',maxIter,...
            'showTrigger',showTrigger,'theta',1,'overRelax',overRelax);
  
  muFit_mVer = muFit2D_mVer( I, z, z0, zR );
  muFit_gBlur = muFit2D_mVer_gBlur( I, z, z0, zR );
  muFit_TV = muFit2D_TV( I, z, z0, zR );
  %muFit_RICR = muFitCVX( lineI, lineMask, z, z0, zR, eta );
  muFit_whTV = muFit2D_whTV( I, z, z0, zR, mask );
  muFit_wTV = weightedTvDenoise_CP( muFit_mVer, paramsCP );

  figure;  hold on;
  %col = round( N / 2 );
  plot( z, trueMu(:,col), 'k', 'LineWidth', 2 );
  plot( z, muFit_mVer(:,col), 'g:', 'LineWidth', 2 );
  plot( z, muFit_gBlur(:,col), 'r--', 'LineWidth', 2 );
  plot( z, muFit_TV(:,col), 'c:', 'LineWidth', 2 );
  plot( z, muFit_whTV(:,col), 'b', 'LineWidth', 2 );
  plot( z, muFit_wTV(:,col), 'k:', 'LineWidth', 2 );
  %plot( z, muFit_RICR, 'c--', 'LineWidth', 2 );
  axis([ 0 max(z) 0 6 ] );
  leg = legend( 'Truth', 'Mod Ver', 'Gauss Smooth', 'TV', 'whTV', 'wTV', ...
    'FontSize', 14, 'Location', 'NorthWest' );
  legH = findobj(leg,'type','text');
  set(legH,'FontSize',20);
  xlabel( 'Depth (mm)', 'FontSize', 16, 'FontWeight', 'bold');
  ylabel( 'Attenuation Coefficient (mm^{-1})', 'FontSize', 16, 'FontWeight', 'bold');


  figure;  hold on;
  I_true = mu2I( trueMu(:,col), z, z0, zR, muAlpha, muBeta, muL0 );
  plot( z, I_true, 'k', 'LineWidth', 2 );
  I_mVer = mu2I( muFit_mVer(:,col), z, z0, zR, muAlpha, muBeta, muL0 );
  plot( z, I_mVer, 'g:', 'LineWidth', 2 );
  I_gBlur = mu2I( muFit_gBlur(:,col), z, z0, zR, muAlpha, muBeta, muL0 );
  plot( z, I_gBlur, 'r--', 'LineWidth', 2 );
  I_whTV = mu2I( muFit_whTV(:,col), z, z0, zR, muAlpha, muBeta, muL0 );
  plot( z, I_whTV, 'b', 'LineWidth', 2 );
  xlabel( 'Depth (mm)', 'FontSize', 16, 'FontWeight', 'bold' );
  ylabel( 'Intensity', 'FontSize', 16, 'FontWeight', 'bold' );
  
  
  trueMuLine = trueMu(:,col);

  vdOffsetThreshPercent = 0.05;
  vD_mVer = findViolationDepth( vdOffsetThreshPercent, trueMuLine, muFit_mVer(:,col), z );
  disp(['Depth of Violation of Modified Vermeer: ', num2str(vD_mVer), 'mm']);

  vD_gBlur = findViolationDepth( vdOffsetThreshPercent, trueMuLine, muFit_gBlur(:,col), z );
  disp(['Depth of Violation of Gaussian Blur: ', num2str(vD_gBlur), 'mm']);

  vD_TV = findViolationDepth( vdOffsetThreshPercent, trueMuLine, muFit_TV(:,col), z );
  disp(['Depth of Violation of Total Variation: ', num2str(vD_TV), 'mm']);

  vD_whTV = findViolationDepth( vdOffsetThreshPercent, trueMuLine, muFit_whTV(:,col), z );
  disp(['Depth of Violation of Weighted Horizontal Total Variation: ', num2str(vD_whTV), 'mm']);

  vD_wTV = findViolationDepth( vdOffsetThreshPercent, trueMuLine, muFit_wTV(:,col), z );
  disp(['Depth of Violation of Weighted Total Variation: ', num2str(vD_wTV), 'mm']);

  etbOffsetThreshPercent = 0.02;
  [etb_mVer, metric_mVer] = findErrorEnergyTooBigDepth( etbOffsetThreshPercent, trueMuLine, muFit_mVer(:,col), z );
  disp(['Depth of Energy Too Big of Modified Vermeer: ', num2str(etb_mVer), 'mm']);

  [etb_gBlur, metric_gBlur] = findErrorEnergyTooBigDepth( etbOffsetThreshPercent, trueMuLine, muFit_gBlur(:,col), z );
  disp(['Depth of Energy Too Big of Gaussian Blur: ', num2str(etb_gBlur), 'mm']);

  [etb_TV, metric_TV] = findErrorEnergyTooBigDepth( etbOffsetThreshPercent, trueMuLine, muFit_TV(:,col), z );
  disp(['Depth of Energy Too Big of TV: ', num2str(etb_TV), 'mm']);

  [etb_whTV, metric_whTV] = findErrorEnergyTooBigDepth( etbOffsetThreshPercent, trueMuLine, muFit_whTV(:,col), z );
  disp(['Depth of Energy Too Big of whTV: ', num2str(etb_whTV), 'mm']);

  [etb_wTV, metric_wTV] = findErrorEnergyTooBigDepth( etbOffsetThreshPercent, trueMuLine, muFit_wTV(:,col), z );
  disp(['Depth of Energy Too Big of wTV: ', num2str(etb_wTV), 'mm']);
end


function [etbDepth, metric] = findErrorEnergyTooBigDepth( offsetThreshPercent, trueMu, muFit, z )
  diff = abs( muFit - trueMu );
  M = numel(muFit);
  errEnergy = zeros( M, 1 );
  for m=1:M
    errEnergy(m) = norm( diff(1:m), 2 );
  end
  metric = errEnergy ./ norm( trueMu, 2 );
  
  etbIndx = find( metric > offsetThreshPercent, 1, 'first' );
  etbDepth = z( etbIndx );
end


function violationDepth = findViolationDepth( offsetThreshPercent, trueMu, muFit, z )
  offsetThresh = trueMu * offsetThreshPercent;
  diff = abs( muFit - trueMu );
  %Offset by 10 to account for Gaussian Blur Kernel
  violations = diff > offsetThresh;
  violationIndx = find( violations(10:end), 1, 'first' );
  violationIndx = violationIndx + 9;
  violationDepth = z( violationIndx );
end

