
function [mu, fos, relFos] = muFit1D_CP2( I, mask, z, z0, zR, eta, muStar )

  fos = [];
  relFos = [];

  nIter = 2000;

  dz = z(2) - z(1);

  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );

  % Determine the norms of A and D
  A = makeA_1D(I, dz, g );
  %normA = powerIterateA1D( rand(size(I)), I, mask, z, dz, z0, zR );
  normA = norm(A);
  Ahat = A / normA;
  normD = 2 / dz;
  normKhat = 1;

  % Choose sigma and tau
  tau = 1 / normKhat;
  sigma = 1 / normKhat;
  theta = 1;

  %Initialize gamma, q, and gammaBar
  gamma = zeros( size(I) );
  qA = zeros( size(I) );
  qDz = zeros( size(I) );
  gammaBar = zeros(size(I));

  b = makeb(I);
  bHat = b / normA;

  M = numel(I);
  Dz = makeD_1D( M, ones(M,1), dz );
  Dhat = Dz / normD;
  fos = zeros(nIter,1);
  
  diffs = zeros(nIter,1);
  for n=1:nIter
    if mod(n,100) == 0
      disp(['Working on CP Iteration: ', num2str(n), ' - diff: ', num2str(diffs(n-1))]);
    end

    fos(n) = objFunction(gamma, I, dz, g, Dz, eta);

    AhatGammaBar = Ahat * gammaBar;
    DhatGammaBar = Dhat * gammaBar;

    %AhatGammaBar = applyAhat( gammaBar, mask, I, dz, z, z0, zR, normA );
    %DhatGammaBar = applyDhat( gammaBar, mask, dz, 1, normD );

    [qA, qDz] = proxFhatConj1D( qA+sigma*AhatGammaBar, ...
      qDz+sigma*DhatGammaBar, bHat, sigma, eta, normA, normD );

    adjAhatQ = Ahat' * qA;
    adjDhatQ = Dhat' * qDz;
    adjKhatQ = adjAhatQ + adjDhatQ;

    oldGamma = gamma;

    gamma = proxG( oldGamma - tau*adjKhatQ );

    gammaBar = gamma + theta*(gamma - oldGamma);

    diffs(n) = norm( gamma - oldGamma );
  end

  mu = 1 ./ gamma;
  muZeroIndxs = find( mask==0 );
  mu( muZeroIndxs ) = 0;

  % convert fos to relative errors
  fStar = objFunction(1./muStar, I, dz, g, Dz, eta);
  relFos = ( fos - fStar ) / fStar;

end
