function [gamma, nIter, relError] = conjGrad_2D( gammaGuess, ...
  I, dz, dx, g, u, yA, yDz, yDx, sysRho, ...
  lambda1, lambda2A, lambda2Dz, lambda2Dx)

  tolerance = 1d-3;
  maxIter = 1000;

  adjKy = applyAdjointK(yA, yDz, yDx, I, dz, dx, g);
  adjKlambda2 = applyAdjointK(lambda2A, lambda2Dz, lambda2Dx,I, dz, dx, g);

  b = u + adjKy - (lambda1/sysRho) - (1/sysRho)*(adjKlambda2);
  normb = norm(b,2);

  %Initialize variables
  warmStart = 1;
  if warmStart==0
    %r = b;
    %gamma = zeros(numel(r), 1);
  else
    gamma = gammaGuess;
    [KgammaA, KgammaDz, KgammaDx] = applyK(gamma, I, dz, dx, g);
    KtKgamma = applyAdjointK(KgammaA, KgammaDz, KgammaDx, I, dz, dx, g);
    r = b - gamma - KtKgamma;
  end
  rho = [Inf sum(r(:)'*r(:)) 0];  %rho[k-2 k-1 k]

  relErrors = zeros(maxIter,1);
  relError = tolerance + 1;
  iter = 0;

  while relError > tolerance && iter < maxIter
    iter = iter + 1;
    if(iter == 1)
      p = r;
    else
      p = r + (rho(2)/rho(1))*p;  % p = r + (rho_{k-1} / rho_{k-2})*p
    end

    [KpA, KpDz, KpDx] = applyK(p, I, dz, dx, g);
    KtKp = applyAdjointK(KpA, KpDz, KpDx, I, dz, dx, g);

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
