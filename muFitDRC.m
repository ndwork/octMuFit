
function mu = muFitDRC( I_in, z_in, h_in )

  goodIndx = find( I_in ~= 0, 1, 'first' );
  I = I_in(goodIndx:end);
  z = z_in(goodIndx:end);
  h = h_in(goodIndx:end);

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

  %x = A \ b;
  pinvA = pinv(A);
  x = pinvA * b;
  muInv = x(2:end);
  mu = zeros( size( I_in ) );
  mu(goodIndx:end) = 1 ./ muInv;
  %mu = 1 ./ x(2:end);
  mu( ~isfinite(mu) ) = 0;
end

