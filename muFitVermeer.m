
function mu = muFitVermeer( I, z )

  mu = zeros( numel(I), 1 );

  dz = z(2)-z(1);

  %Itilde = I ./ g;
  Itilde = I;
  
  n = numel(I);
  ut = triu( ones(n) );
  
  mu = 1/(2*dz) * log( 1 + Itilde ./ (ut*I) );

end


