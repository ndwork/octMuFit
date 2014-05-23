function [gamma] = conjGrad_2D_MATLAB( gammaGuess, I, z, dz, dx, z0, zR, u, yA, yDz, yDx, sysRho, lambda1, lambda2A, lambda2Dz, lambda2Dx)
  
  % Want to solve (I + K'K)gamma = u + K'y - lambda1/rho - K'*lambda2/rho
  % for gamma
  
  adjKy = applyAdjointK(yA, yDz, yDx, I, z, dz, dx, z0, zR);
  adjKlambda2 = applyAdjointK(lambda2A, lambda2Dz, lambda2Dx,I, z, dz, dx, z0, zR);
  
  b = u + adjKy - (lambda1/sysRho) - (1/sysRho)*(adjKlambda2);
      
  eps = 0.5;
  nIter = 100;
  
  function [Ap] = applyA_2DconjGrad(p)  
    % Function to calculate (I + K'K)*p
    % Function is nested so that it's scope includes the variables passed
    % in to conjGrad_2D_MATLAB
    
    % Need to reshape p so that's m by n (size of I)
    p = reshape(p, size(I));
    
    [KpA, KpDz, KpDx] = applyK(p, I, z, dz, dx, z0, zR);
    KtKp = applyAdjointK(KpA, KpDz, KpDx, I, z, dz, dx, z0, zR);
    Ap = p + KtKp;   
    
    % Need to reshape Ap so that its mn by 1 because cgs only works with 1D
    % problems
    Ap = reshape(Ap, numel(p), 1);
  end
  
  [gamma,flag,relres,iter] = cgs(@applyA_2DconjGrad, b, 1d-3, nIter, [], [], reshape(gammaGuess, [numel(I), 1]));
  
  
  
end
