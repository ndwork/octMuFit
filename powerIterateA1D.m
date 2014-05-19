

function normA = powerIterateA1D( x, I, mask, z, dz, z0, zR )
  % Use power iteration to estimate the norm of K
  
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  A = makeA_1D( I, dz, g );
  %K = makeK_1D( I, dz, g, mask );

  tolerance = 1d-5;
  diff = tolerance+1;

  maxNIter = 1000;

  lambda = 0;
  iter = 0;
  while diff > tolerance && iter < maxNIter
    iter = iter + 1;
    disp(['Working on A Power Iteration: ', num2str(iter)]);

    %yA = applyA( gamma, mask, I, dz, z, z0, zR );
    %AtAx = applyAdjointA( yA, mask, I, z, dz, z0, zR );
    
    AtAx = A' * A * x;

    lastLambda = lambda;
    lambda = norm( AtAx(:), 2 );
    x = AtAx / lambda;

    diff = abs( lambda - lastLambda );
  end

  normA = sqrt(lambda);

end

