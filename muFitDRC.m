
function mu = muFitDRC( I, z, h )

  dz = z(2) - z(1);
  Ihat = I ./ h;

  %f = 1 ./ h .* ( I.*I ./ (I.*I + (1d-3)^2 ) );
  %Ihat = I .* f;

  b = cumsum( Ihat ) * dz;

  M = numel(z);
  A = zeros(M+1);

  A(1:M,1) = 1;
  for i=1:M
    A(i,i+1) = -0.5 * Ihat(i);
  end

  % Additional constraint
  b = [ b; 0 ];
  A(end,:) = [ zeros(1,M-1) 1 -1 ];

  x = A \ b;
  muInv = x(2:end);
  mu = 1 ./ muInv;
  mu( ~isfinite(mu) ) = 0;
end

