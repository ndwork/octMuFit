
function h = makeConfocalFunction( z, z0, zR )
  % z0 is the focal plane location

  tmp = (z-z0)/zR;
  h = 1 ./ ( tmp.^2 + 1 );

[lambda,deltaLambda,dLambda] = getTelestoFalloffParams();
f = makeFalloffFunction( z, lambda, deltaLambda, dLambda );
h = h .* f;

end
