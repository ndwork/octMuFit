
function mu = muFitModVermeer( I, z, h )

  dz = z(2) - z(1);
  term = I./h;
  integral = zeros(size(term));
  nZ = numel(z);
  for m = 1:(nZ - 1)
     integral(m) = sum(term(m:end - 1)); 
  end

  mu = (.5/dz)*(term./integral);
  %mu = 0.5/dz * log( 1 + term./integral );

end

