
function normK = powerIterateK1D( x, I, z, dz, z0, zR )
  % Use power iteration to estimate the norm of K
  
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  %K = makeK_1D( I, dz, g, mask );
  
  tolerance = 1d-5;
  diff = tolerance+1;
  
  maxNIter = 1000;
  
  lambda = 0;
  iter = 0;
  while diff > tolerance && iter < maxNIter
    iter = iter + 1;
    disp(['Working on Power Iteration: ', num2str(iter)]);

    [yA, yDz] = applyK1D( x, I, z, dz, z0, zR );
    KtKx = applyAdjointK1D( yA, yDz, I, z, dz, z0, zR );

    %KtKx = K' * K * x;

    lastLambda = lambda;
    lambda = norm( KtKx(:), 2 );
    x = KtKx / lambda;

    diff = abs( lambda - lastLambda );
  end

  normK = sqrt(lambda);

end

