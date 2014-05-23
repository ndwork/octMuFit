function [gamma] = conjGrad_1D_MATLAB( gammaGuess, I, z, dz, z0, zR, u, yA, yDz, sysRho, lambda1, lambda2A, lambda2Dz)
  
  % Want to solve (I + K'K)gamma = u + K'y - lambda1/rho - K'*lambda2/rho
  % for gamma where I and gamma are 1D
  
   adjKy = applyAdjointK1D(yA, yDz, I, z, dz, z0, zR);
   adjKlambda2 = applyAdjointK1D(lambda2A, lambda2Dz, I, z, dz, z0, zR);
   
   b = u + adjKy - (lambda1/sysRho) - (1/sysRho)*(adjKlambda2);
   
   eps = 0.5;
   nIter = 100;
   
  function [Ap] = applyA_1DconjGrad(p)
    % Function to calculate (I + K'K)*p
    % Function is nested so that it's scope includes the variables passed
    % in to conjGrad_1D_MATLAB
    [KpA, KpDz] = applyK1D(p, I, z, dz, z0, zR);
    KtKp = applyAdjointK1D(KpA, KpDz, I, z, dz, z0, zR);
    Ap = p + KtKp;
  end
  
  [gamma,flag,relres,iter] = cgs(@applyA_1DconjGrad, b, 1d-3, nIter, [], [], gammaGuess);  
   
  
end
