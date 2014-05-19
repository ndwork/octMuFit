function [mu, fos, relFos] = muFit1D_ADMMwCVX( I, mask, z, z0, zR, eta, muStar )
  
  nIter = 10;

  dz = z(2)-z(1);

  %Choose rho
  rho = 1d4;

  %Initialize gamma, q, and u
  gamma = zeros( size(I) );
  zA = zeros( size(I) );
  zDz = zeros( size(I) );
  
  
  yA = zeros( size(I) );
  yDz = zeros( size(I) );
  y = [yA; yDz];
  
  b = makeb(I);
  
  % Make g and D for objective eval
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  Dz = makeD_1D(length(gamma), mask, dz);
  fos = zeros(nIter, 1);
  
  K = makeK_1D(I, dz, g, mask);

  diffs = zeros(nIter, 1);
  for n = 1:nIter
    if mod(n,2)==0
      disp(['Working on 1D ADMM Iteration: ', num2str(n), ' - diff: ', num2str(diffs(n-1))]);
    end
    
    fos(n) = objFunction(gamma, I, dz, g, Dz, eta);   
    
    oldGamma = gamma;
    
    cvx_begin quiet
      variable CVX_gamma(size(I))
      minimize (y'*(K*CVX_gamma - [zA; zDz]) + (rho/2)*sum((K*CVX_gamma - [zA; zDz]).^2))
      CVX_gamma >= 0
    cvx_end
    
    gamma = CVX_gamma;
    
    cvx_begin quiet
      variables CVX_zA(size(I)) CVX_zDz(size(I))
      minimize (0.5*sum((CVX_zA - b).^2) + eta*sum(abs(CVX_zDz)) + y'*(K*CVX_gamma - [zA; zDz]) + (rho/2)*sum((K*CVX_gamma - [zA; zDz]).^2))
    cvx_end
    
    zA = CVX_zA;
    zDz = CVX_zDz;
    
    y = y + rho*(K*gamma - [zA; zDz]);   
    
    diffs(n) = norm( gamma - oldGamma );
  end
  mu = 1 ./ gamma;
  muZeroIndxs = find( mask==0 );
  mu( muZeroIndxs ) = 0;

  % convert fos to relative errors
  fStar = objFunction(1./muStar, I, dz, g, Dz, eta);
  relFos = ( fos - fStar ) / fStar;
end