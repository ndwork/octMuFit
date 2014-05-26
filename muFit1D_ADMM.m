
function [mu, fos, relFos] = muFit1D_ADMM(I, mask, z, z0, zR, eta, muStar )

  fos = [];
  relFos = [];

  rho = 1d-3;
  nIter = 400;

  M = numel(I);
  dz = z(2) - z(1);

  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );

  %A = makeA_1D( I, dz, g );
  D = makeD_1D( M, dz );
  b = makeb( I );

  u = zeros(M,1);
  y = zeros(2*M,1);
  gamma = zeros(M,1);
  lambda1 = zeros(M,1);
  lambda2 = zeros(2*M,1);

  fos = zeros(nIter,1);

  %K = [A; D];
  %IKtK = eye(M) + K'*K;

  for n=1:nIter
    if mod(n,50)==0 disp(['1D ADMM iteration: ', num2str(n)]); end;
    fos(n) = objFunction(gamma, I, mask, dz, g, D, eta);

    %gamma = IKtK \ ( u + K'*y - 1/rho*lambda1 - 1/rho*K'*lambda2 );
    if( n > 1 )
      gamma = conjGrad_1D(gamma, I, dz, g, u, y(1:M), ...
        y(M+1:end), rho, lambda1, lambda2(1:M), lambda2(M+1:end));
    end

    %Agamma = A*gamma;
    %Dgamma = D*gamma;
    Agamma = applyA( gamma, I, dz, g );
    Dgamma = applyD( gamma, dz, 1 );
    
    [y1, y2] = proxF1D( Agamma + 1/rho*lambda2(1:M), ...
                        Dgamma + 1/rho*lambda2(M+1:end), ...
                        mask, b, 1/rho, eta);
    %y = [ y1; y2 ];

    u = proxG(gamma + 1/rho*lambda1, mask);

    lambda1 = lambda1 + rho*(gamma - u);
    lambda2(1:M) = lambda2(1:M) + rho*(Agamma-y1);
    lambda2(M+1:end) = lambda2(M+1:end) + rho*(Dgamma-y2);
  end

  mu = 1 ./ gamma;
  mu( mask==0 ) = 0;

  % convert fos to relative errors
  fStar = objFunction(1./muStar, I, mask, dz, g, D, eta);
  relFos = ( fos - fStar ) / fStar;

end


