
function [mu, fos] = muFit2D_ADMM(I, mask, z, dx, z0, zR, etaz, etax )

  fos = [];

  rho = 0.001;
  nIter = 400;

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

  fos = zeros(nIter, 1);
  cgNIters = zeros( nIter, 1 );
  cgRelErrors = zeros( nIter, 1 );
  
  for n=1:nIter
    if mod(n,10)==0 disp(['2D ADMM iteration: ', num2str(n)]); end;

    fos(n) = objFunction2D(gamma, I, mask, dz, dx, z, z0, zR, etaz, etax);

    if( n > 1 )
      [gamma, cgNIter, cgRelError] = conjGrad_2D(gamma, ...
        I, z, dz, dx, z0, zR, u, y1, y2, y3, ...
        rho, lambda1, lambda2(1:M, :), ...
        lambda2(M+1:2*M, :), lambda2(2*M+1 : end, :));
      cgNIters(n) = cgNIter;
      cgRelErrors(n) = cgRelError;
    end

    Agamma = applyA( gamma, I, dz, z, z0, zR );
    DzGamma = applyD( gamma, dz, 1 );
    DxGamma = applyD( gamma, dx, 2 );

    [y1, y2, y3] = proxF2D( Agamma + 1/rho*lambda2(1:M,:), ...
                       DzGamma + 1/rho*lambda2(M+1:2*M,:), ...
                       DxGamma + 1/rho*lambda2(2*M+1:end,:), ...
                       mask, b, 1/rho, etaz, etax );

    u = proxG(gamma + 1/rho*lambda1, mask);

    lambda1 = lambda1 + rho*(gamma - u);
    lambda2(1:M,:) = lambda2(1:M,:) + rho*(Agamma-y1);
    lambda2(M+1:2*M,:) = lambda2(M+1:2*M,:) + rho*(DzGamma-y2);
    lambda2(2*M+1:end,:) = lambda2(2*M+1:end,:) + rho*(DxGamma-y3);
  end

  mu = 1 ./ gamma;
  mu( mask==0 ) = 0;
  
  figure, plot( cgNIters );  title('CG N Iterations');
  figure, plot( cgRelErrors );  title('CG Relative Errors');
  
end
