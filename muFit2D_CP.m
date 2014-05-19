function mu = muFit2D_CP( I, mask, z, dz, dx, z0, zR )
  % Use Chambolle Pock to extract mu form 2D OCT data

  I = I ./ max(I(:));
  
  nIter = 10000;

  % Determine the norm of K
  normK = powerIterateK( rand(size(I)), I, mask, z, dz, dx, z0, zR );

  % Choose sigma and tau
  tau = 1/sqrt(normK);
  sigma = 1/(tau*normK^2);

  %Initialize gamma, q, and gammaBar
  gamma = zeros( size(I) );
  qA = zeros( size(I) );
  qDx = zeros( size(I) );
  qDz = zeros( size(I) );

  gammaBar = zeros(size(I));

  b = makeb(I);

  eta = 100;
  theta = 1;
  diffs = zeros(nIter,1);
  for n=1:nIter
    if mod(n,10)==0
      disp(['Working on CP Iteration: ', num2str(n), ' - diff: ', num2str(diffs(n-1))]);
    end

    [AgammaBar, DzGammaBar, DxGammaBar] = applyK( gammaBar, mask, I, z, dz, dx, z0, zR );

    [qA qDz qDx] = proxFconj(qA+sigma*AgammaBar, qDz+sigma*DzGammaBar, qDx + sigma*DxGammaBar, b, sigma, eta);

    adjKGamma  = applyAdjointK(qA, qDz, qDx, I, mask, z, dz, dx, z0, zR);

    oldGamma = gamma;

    gamma = proxG(oldGamma - tau*adjKGamma);

    gammaBar = gamma + theta*(gamma - oldGamma);

    diffs(n) = norm( gamma - oldGamma );
  end

  mu = 1 ./ gamma;
  muZeroIndxs = find( mask==0 );
  mu( muZeroIndxs ) = 0;

end

