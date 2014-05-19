
function [mu, fos, relFos] = muFit1D_CP( I, mask, z, z0, zR, eta, muStar, paramsCP )

  nIter = paramsCP.nIter;
  tau = paramsCP.tau;  
  
  dz = z(2) - z(1);

  % Determine the norm of K
%   normK = powerIterateK1D( rand(size(I)), I, mask, z, dz, z0, zR ); 

  % Choose sigma and tau
%   tau = 1 / normK;
%   sigma = 1 / normK;
    
%   sigma = 1/(tau*normK^2);
  theta = 1;

  %Initialize gamma, q, and gammaBar
  gamma = zeros( size(I) );
  qA = zeros( size(I) );
  qDz = zeros( size(I) );
  gammaBar = zeros(size(I));

%   b = makeb(I);
  
  %%%
  b = 1:numel(I);
  b = b'/10;
  b = sin(b);
  b = b + .05*randn(numel(I),1);
  
  %%%

  M = numel(I);
  Dz = makeD_1D( M, ones(M,1), dz );
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  fos = zeros(nIter,1);
  
%   appA = @(u) applyA( u, mask, I, dz, z, z0, zR );
%   appD = @(u) applyD( u, mask, dz, 1 );  
  
  appD = @(u) applyD(u,mask,1,1);
  
  appA = @(u) u; 
  appK = @(u) [appA(u);appD(u)];
  
%   appATrans = @(v) applyAdjointA( v , mask, I, z, dz, z0, zR );
%   appDTrans = @(v) applyAdjointD( v, mask, dz, 1 );
  appATrans = @(v) v;
  appDTrans = @(v) applyAdjointD( v, mask, 1, 1 );
  appKTrans = @(w) appATrans(w(1:M)) + appDTrans(w(M+1:end));   
  
  nrmK = estimateNormByPowerIteration(appK,appKTrans,randn(M,1),1000);
  sigma = 1/(tau*nrmK^2);

  diffs = zeros(nIter,1);
  for n=1:nIter
    if mod(n,100) == 0
      disp(['Working on CP Iteration: ', num2str(n), ' - diff: ', num2str(diffs(n-1))]);      
    end

%     fos(n) = objFunction(gamma, I, dz, g, Dz, eta);
    fos(n) = .5*sum( (appA(gamma) - b).^2 ) + eta*sum(abs(appD(gamma)));    

%     [AgammaBarCheck, DzGammaBarCheck] = applyK1D( gammaBar, mask, I, z, dz, z0, zR );
    AgammaBar = appA(gammaBar);
    DzGammaBar = appD(gammaBar);
    
    z1 = qA + sigma*AgammaBar; 
    z2 = qDz + sigma*DzGammaBar;
        
%     [qACheck, qDzCheck] = proxFconj1D(qA+sigma*AgammaBar, qDz+sigma*DzGammaBar, b, sigma, eta);
        
    qA = (z1 - sigma*b)/(1 + sigma);     
        
    qDz = min(z2/eta,1);
    qDz = max(qDz,-1);
    qDz = eta*qDz;            

%     adjKqCheck  = applyAdjointK1D(qA, qDz, I, mask, z, dz, z0, zR);    
    
    adjKq = appKTrans([qA;qDz]);

    oldGamma = gamma;

%     gammaCheck = proxG( oldGamma - tau*adjKq );
    gamma = max(oldGamma - tau*adjKq, 0);    

    gammaBar = gamma + theta*(gamma - oldGamma);

    diffs(n) = norm( gamma - oldGamma );
  end

% plot( b );
% hold on;
% A = makeA_1D( I, dz, g );
% plot( A*gamma, 'r' );
% ylim([0 max(b)]);
% disp(norm(A*gamma-b));

  
  mu = 1 ./ gamma;
  muZeroIndxs = find( mask==0 );
  mu( muZeroIndxs ) = 0;

  % convert fos to relative errors
  fStar = objFunction(1./muStar, I, dz, g, Dz, eta);
  relFos = ( fos - fStar ) / fStar;

end





