
function [mu, fos] = muFit2D_ADMM(I, mask, z, dx, z0, zR, etaz, etax )

  fos = [];

  rho = 0.001;
  nIter = 800;

  [M, N] = size(I);
  dz = z(2) - z(1);


  %normK = powerIterateK( rand(size(I)), I, z, dz, dx, z0, zR );

  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );

  b = makeb( I );

  u = zeros(M,N);
  gamma = zeros(M,N);
  lambda1 = zeros(M,N);
  lambda2 = zeros(3*M,N);
  
  uOld = zeros(size(u));  %Initialize for residual calculation
  y1Old = zeros(M, N);
  y2Old = zeros(M, N);
  y3Old = zeros(M, N);


  fos = zeros(nIter, 1);
  cgNIters = zeros( nIter, 1 );
  cgRelErrors = zeros( nIter, 1 );
  
  epsRel = 1d-3;
  epsAbs = 0;
  
  primalRes = zeros(nIter, 1);
  dualRes = zeros(nIter, 1);
  epsPrimal = zeros(nIter, 1);
  epsDual = zeros(nIter, 1);
  
  for n=1:nIter
    if mod(n,10)==0 disp(['2D ADMM iteration: ', num2str(n)]); end;

    fos(n) = objFunction2D(gamma, I, mask, dz, dx, g, etaz, etax);

    if( n > 1 )
      [gamma, cgNIter, cgRelError] = conjGrad_2D(gamma, ...
        I, dz, dx, g, u, y1, y2, y3, ...
        rho, lambda1, lambda2(1:M, :), ...
        lambda2(M+1:2*M, :), lambda2(2*M+1 : end, :));
      cgNIters(n) = cgNIter;
      cgRelErrors(n) = cgRelError;
    end

    Agamma = applyA( gamma, I, dz, g );
    DzGamma = applyD( gamma, dz, 1 );
    DxGamma = applyD( gamma, dx, 2 );
    
    if(n > 1)
      y1Old = y1;
      y2Old = y2;
      y3Old = y3;
    end
    [y1, y2, y3] = proxF2D( Agamma + 1/rho*lambda2(1:M,:), ...
                       DzGamma + 1/rho*lambda2(M+1:2*M,:), ...
                       DxGamma + 1/rho*lambda2(2*M+1:end,:), ...
                       mask, b, 1/rho, etaz, etax );

    uOld = u;
    u = proxG(gamma + 1/rho*lambda1, mask);

    lambda1 = lambda1 + rho*(gamma - u);
    lambda2(1:M,:) = lambda2(1:M,:) + rho*(Agamma-y1);
    lambda2(M+1:2*M,:) = lambda2(M+1:2*M,:) + rho*(DzGamma-y2);
    lambda2(2*M+1:end,:) = lambda2(2*M+1:end,:) + rho*(DxGamma-y3);
    
    % Determine primal and dual feasibiilty
    % See Distributed Optimization and Statistical Learning via ADMM pg 19    
    norm_rk = sum(Agamma(:).^2 + DzGamma(:).^2 + DxGamma(:).^2)...
                + sum(u(:).^2 + y1(:).^2 + y2(:).^2 + y3(:).^2)...
                - 2*(gamma(:)'*u(:) + Agamma(:)'*y1(:) + DzGamma(:)'*y2(:) + DxGamma(:)'*y3(:)); 
              
    normAx = norm(gamma(:)) + norm(Agamma(:)) + norm(DzGamma(:)) + norm(DxGamma(:));
    normBz = norm(u(:)) + norm(y1(:)) + norm(y2(:)) + norm(y3(:));
    normc = 0;
    epsPrimal(n) = (sqrt(M*N)*epsAbs + epsRel *max([normAx, normBz, normc])).^2;
    
    yDiff1 = y1 - y1Old;
    yDiff2 = y2 - y2Old;
    yDiff3 = y3 - y3Old;
    uDiff = u - uOld;
    
    adjKyDiff = applyAdjointK(yDiff1, yDiff2, yDiff3, I, dz, dx, g);
    adjKlambda2 = applyAdjointK(lambda2(1:M,:), lambda2(M+1:2*M,:), lambda2(2*M+1:end,:), I, dz, dx, g);
    
    norm_sk = norm(rho*adjKyDiff(:)) + norm(rho*uDiff(:));    
    
    epsDual(n) = (sqrt(4*M)*epsAbs + epsRel*(norm(lambda1(:)) + norm(adjKlambda2(:)))).^2;
    
    primalRes(n) = norm_rk;
    dualRes(n) = norm_sk;
  end

  mu = 1 ./ gamma;
  mu( mask==0 ) = 0;
  
  figure, plot( cgNIters );  title('CG N Iterations');
  figure, plot( cgRelErrors );  title('CG Relative Errors');
  figure, semilogy(primalRes), hold all, semilogy(epsPrimal);
  legend('||r^k||_2^2','\epsilon_{primal}')
  figure, semilogy(dualRes), hold all, semilogy(epsDual);
  legend('||s^k||_2^2', '\epsilon_{dual}')

end
