
function normK = powerIterateK( x, I, mask, z, dz, dx, z0, zR )
  % Use power iteration to estimate the norm of K
  
  tolerance = 1d-5;
  diff = tolerance+1;
  
  maxNIter = 1000;
  
  lambda = 0;
  iter = 0;
  while diff > tolerance && iter < maxNIter
    iter = iter + 1;
    disp(['Working on Power Iteration: ', num2str(iter)]);

    [yA, yDz, yDx] = applyK( x, mask, I, z, dz, dx, z0, zR );
    KtKx = applyAdjointK( yA, yDz, yDx, I, mask, z, dz, dx, z0, zR );

    lastLambda = lambda;
    lambda = norm( KtKx(:), 2 );
    x = KtKx / lambda;

    diff = abs( lambda - lastLambda );
  end

  normK = sqrt(lambda);

end

