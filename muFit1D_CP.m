
function [mu, fos, relFos] = muFit1D_CP( I, mask, z, z0, zR, eta, muStar )

  nIter = 2000;

  dz = z(2) - z(1);

  % Determine the norm of K
  normK = powerIterateK1D( rand(size(I)), I, mask, z, dz, z0, zR );

  % Choose sigma and tau
  tau = 1 / normK;
  sigma = 1 / normK;
  theta = 1;

  %Initialize gamma, q, and gammaBar
  gamma = zeros( size(I) );
  qA = zeros( size(I) );
  qDz = zeros( size(I) );
  gammaBar = zeros(size(I));

  b = makeb(I);

  M = numel(I);
  Dz = makeD_1D( M, ones(M,1), dz );
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  fos = zeros(nIter,1);

  diffs = zeros(nIter,1);
  for n=1:nIter
    if mod(n,100) == 0
      disp(['Working on CP Iteration: ', num2str(n), ' - diff: ', num2str(diffs(n-1))]);
    end

    fos(n) = objFunction(gamma, I, dz, g, Dz, eta);

    [AgammaBar, DzGammaBar] = applyK1D( gammaBar, mask, I, z, dz, z0, zR );

    [qA, qDz] = proxFconj1D(qA+sigma*AgammaBar, qDz+sigma*DzGammaBar, b, sigma, eta);

    adjKq  = applyAdjointK1D(qA, qDz, I, mask, z, dz, z0, zR);

    oldGamma = gamma;

    gamma = proxG( oldGamma - tau*adjKq );

    gammaBar = gamma + theta*(gamma - oldGamma);

    diffs(n) = norm( gamma - oldGamma );
  end

% plot( b );
% hold on;
% A = makeA_1D( I, dz, g );
% plot( A*gamma, 'r' );
% ylim([0 max(b)]);
% disp(norm(A*gamma-b));

  
  mu = 1 ./ gamma;
  muZeroIndxs = find( mask==0 );
  mu( muZeroIndxs ) = 0;

  % convert fos to relative errors
  fStar = objFunction(1./muStar, I, dz, g, Dz, eta);
  relFos = ( fos - fStar ) / fStar;

end





