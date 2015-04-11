
function mu = muFitModVermeer( I, z, h )

  dz = z(2) - z(1);

  %term = I ./ h;

  f = 1 ./ h .* ( I.*I ./ (I.*I + (0.5d-2)^2 ) );  % for everything else
  %f = 1 ./ h .* ( I.*I ./ (I.*I + (1d-1)^2 ) );  % for intralipid
  term = I .* f;

  nTerm = numel(term);
  integral = zeros(nTerm,1);
  for m = 1:(nTerm - 1)
    integral(m) = sum( term(m:end-1) );
  end

  mu = (.5/dz)*(term./integral);
  %mu = 0.5/dz * log( 1 + term./integral );
  
  mu( ~isfinite(mu) ) = 0;
end

