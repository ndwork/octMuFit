function mu = muFit2D_LADMM( I, mask, z, dz, dx, z0, zR )
  % Use Linear ADMM to extract mu from 2D OCT data

  nIter = 1000;

  % Determine the norm of K
  normK = powerIterateK( rand(size(I)), I, mask, z, dz, dx, z0, zR );

  % Choose sigma and tau
  tau = 1;
  lambda = normK*tau;

  %Initialize gamma, q, and u
  gamma = zeros( size(I) );
  qA = zeros( size(I) );
  qDx = zeros( size(I) );
  qDz = zeros( size(I) );

  uA = zeros( size(I) );
  uDx = zeros( size(I) );
  uDz = zeros( size(I) );

  b = makeb(I);

  eta = 100;
  diffs = zeros(nIter, 1);
  for n = 1:nIter
    if n > 1
      disp(['Working on ADMM Iteration: ', num2str(n), ' - diff: ', num2str(diffs(n-1))]);
    end
    
    [Agamma, DzGamma, DxGamma] = applyK( gamma, mask, I, z, dz, dx, z0, zR );
    
    tmpA = Agamma - qA + uA;
    tmpDz = DzGamma - qDz + uDz;
    tmpDx = DxGamma - qDx + uDx;
    
    adjKtmp = applyAdjointK( tmpA, tmpDz, tmpDx, I, mask, z, dz, dx, z0, zR );
    
    oldGamma = gamma;
    
    gamma = proxG(gamma - (tau/lambda) * adjKtmp);
    
    [Agamma, DzGamma, DxGamma] = applyK( gamma, mask, I, z, dz, dx, z0, zR );
    
    [qA, qDz, qDx] = proxF(Agamma + uA, DzGamma + uDz, DxGamma + uDx, b, lambda, eta);    
    
    uA = uA + Agamma - qA;
    uDz = uDz + DzGamma - qDz;
    uDx = uDx + DxGamma - qDx;
    
    diffs(n) = norm( gamma - oldGamma );
  end
  
  mu = 1 ./ gamma;
  muZeroIndxs = find( mask==0 );
  mu( muZeroIndxs ) = 0;

end