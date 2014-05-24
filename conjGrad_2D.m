function [gamma, nIter, relError] = conjGrad_2D( gammaGuess, ...
  I, z, dz, dx, z0, zR, u, yA, yDz, yDx, sysRho, ...
  lambda1, lambda2A, lambda2Dz, lambda2Dx)

  adjKy = applyAdjointK(yA, yDz, yDx, I, z, dz, dx, z0, zR);
  adjKlambda2 = applyAdjointK(lambda2A, lambda2Dz, lambda2Dx,I, z, dz, dx, z0, zR);

  b = u + adjKy - (lambda1/sysRho) - (1/sysRho)*(adjKlambda2);
  normb = norm(b,2);

  eps = 0.5;
  maxIter = 1000;

  %Initialize variables
  warmStart = 1;
  if warmStart==0
    %r = b;
    %gamma = zeros(numel(r), 1);
  else
    gamma = gammaGuess;
    [KgammaA, KgammaDz, KgammaDx] = applyK(gamma, I, z, dz, dx, z0, zR);
    KtKgamma = applyAdjointK(KgammaA, KgammaDz, KgammaDx, I, z, dz, dx, z0, zR);
    r = b - gamma - KtKgamma;
  end
  rho = [Inf sum(r(:)'*r(:)) 0];  %rho[k-2 k-1 k]

  relErrors = zeros(maxIter,1);
  tolerance = 1d-3;
  relError = tolerance + 1;
  iter = 0;

  while relError > tolerance && iter < maxIter
    iter = iter + 1;
    if(iter == 1)
      p = r;
    else
      p = r + (rho(2)/rho(1))*p;  % p = r + (rho_{k-1} / rho_{k-2})*p
    end

    [KpA, KpDz, KpDx] = applyK(p, I, z, dz, dx, z0, zR);
    KtKp = applyAdjointK(KpA, KpDz, KpDx, I, z, dz, dx, z0, zR);

    w = p + KtKp;    % w = (I + K'*K)*p

    alpha = rho(2)/(p(:)'*w(:)); % alpha = rho_{k-1} / (p'*w)

    gamma = gamma + alpha*p;

    r = r - alpha*w;

    rho(3)= r(:)'*r(:);   % rho_k = r'*r
    rho = circshift(rho, [0, -1] );  % shift the rhos up so that we can update rho_k in the next iteration    

    relError = sqrt( rho(2) ) / normb;
    relErrors(iter) = relError;
  end
  nIter = iter;

end
