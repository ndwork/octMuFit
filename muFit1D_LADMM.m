function [mu, fos, relFos] = muFit1D_LADMM( I, mask, z, z0, zR, eta, muStar )
  % Use Linear ADMM to extract mu from 1D OCT data

  nIter = 10000;

  dz = z(2)-z(1);

  %Initialize gamma, q, and u
  gamma = zeros( size(I) );
  qA = zeros( size(I) );
  qD = zeros( size(I) );
  uA = zeros( size(I) );
  uD = zeros( size(I) );
  
  b = makeb(I);
  
  % Make g and D for objective eval
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  
  fos = zeros(nIter, 1);
  
  D = makeD_1D(length(gamma), mask, dz);
  
  % Choose lambda and tau
  normK = powerIterateK1D(rand(size(gamma), 1), I, mask, z, dz, z0, zR);
  lambda = normK*normK;
  tau = 0.5; 

  diffs = zeros(nIter, 1);
  for n = 1:nIter
    if mod(n,100)==0
      disp(['Working on 1D LADMM Iteration: ', num2str(n), ' - diff: ', num2str(diffs(n-1))]);
    end

    fos(n) = objFunction(gamma, I, dz, g, D, eta);

    [Agamma, Dgamma] = applyK1D(gamma, mask, I, z, dz, z0, zR);
    
    tmpA = Agamma - qA + uA;
    tmpD = Dgamma - qD + uD;

    adjKtmp = applyAdjointK1D(tmpA, tmpD, I, mask, z, dz, z0, zR);

    oldGamma = gamma;
    gamma = proxG(gamma - (tau/lambda) * adjKtmp);

    [Agamma, Dgamma] = applyK1D(gamma, mask, I, z, dz, z0, zR);

    [qA, qD] = proxF1D(Agamma + uA, Dgamma + uD, b, lambda, eta); 

    uA = uA + Agamma - qA;
    uD = uD + Dgamma - qD;

    diffs(n) = norm( gamma - oldGamma );
  end

  mu = 1 ./ gamma;
  muZeroIndxs = find( mask==0 );
  mu( muZeroIndxs ) = 0;

  % convert fos to relative errors
  fStar = objFunction(1./muStar, I, dz, g, D, eta);
  relFos = ( fos - fStar ) / fStar;
end

