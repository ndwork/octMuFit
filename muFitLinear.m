
function mu = muFitLinear( I, z, x0, z0 )
  %  Note that this doesn't enforce the non-negativity constraint

  n = numel(I);
  dz = z(2) - z(1);
  ut = triu( ones(n,n) ) * dz;

  g = 1 ./ sqrt( ( ( z - x0 ) / z0 ).^2 +1 );
  tmp = (z-x0)/z0;
  gPrime = -tmp ./ ( tmp.^2 + 1 ).^(3/2);
  
  A = diag(I)./2 + ut*diag(I./g.*gPrime./2);
  b = ut*I;
  gam = A \ b;
  mu = 1 ./ gam;

  %condNumA = cond( A );
  %disp(['The condition number of A is: ', num2str(condNumA)]);

end
