function [gamma] = conjGrad_1D( gammaGuess, I, z, dz, z0, zR, u, yA, yDz, sysRho, lambda1, lambda2A, lambda2Dz)
  
  adjKy = applyAdjointK1D(yA, yDz, I, z, dz, z0, zR);
  adjKlambda2 = applyAdjointK1D(lambda2A, lambda2Dz, I, z, dz, z0, zR);
  
  b = u + adjKy - (lambda1/sysRho) - (1/sysRho)*(adjKlambda2);
    
  eps = 0.5;
  nIter = 100;

  ident = eye( numel(z) ); 
  
  %Initialize variables
  warmStart = 1;
  if warmStart==0
    %r = b;
    %gamma = zeros(numel(r), 1);
  else
    gamma = gammaGuess;
    [KgammaA, KgammaDz] = applyK1D(gamma, I, z, dz, z0, zR);
    KtKgamma = applyAdjointK1D(KgammaA, KgammaDz, I, z, dz, z0, zR);
    r = b - ident*gamma - KtKgamma;
  end
  rho = [Inf r'*r 0];  %rho[k-2 k-1 k]
    
  for k = 1:nIter
%     if(sqrt(rho(2)) <= eps*norm(b))
%       break;
%     end
    
    if(k == 1)
      p = r;
    else
      p = r + (rho(2)/rho(1))*p;  % p = r + (rho_{k-1} / rho_{k-2})*p
    end
    
    [KpA, KpDz] = applyK1D(p, I, z, dz, z0, zR);
    KtKp = applyAdjointK1D(KpA, KpDz, I, z, dz, z0, zR);
    
    w = ident*p + KtKp;    % w = (I + K'*K)*p
    
    alpha = rho(2)/(p'*w); % alpha = rho_{k-1} / (p'*w)
    
    gamma = gamma + alpha*p;
    
    r = r - alpha*w;
    
    rho(3)= r'*r;   % rho_k = r'*r
    
    rho = circshift(rho, [0, -1] );  % shift the rhos up so that we can update rho_k in the next iteration

    
  end  
  
end