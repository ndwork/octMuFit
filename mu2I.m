
function I = mu2I( mu, z, z0, zR, alpha, beta, L0 )
  % z represents the depth at the front of each dz layer
  % I represents the intensity in the middle of each dz layer

  dz = z(2:end) - z(1:end-1);
  dz = [ dz; dz(end) ];  % symmetric boundary condition

  g = 1 ./ sqrt( ((z + dz/2 - z0) / zR).^2 + 1 );

  I = zeros( numel(z), 1 );

  k = alpha*beta*L0;
  if numel(k) == 0 k=1; end;

  tmp = mu(1) * dz(1)/2;
  for i=1:numel(I)
    I(i) = k*mu(i) .* exp( -2 * tmp ) .* g(i);
    if i < numel(I)
      tmp = tmp + mu(i)*dz(i)/2 + mu(i+1)*dz(i+1)/2;
    end
  end  

%   I(1) = k * mu(1) * g(1);
%   tmp = 0;
%   for i=2:numel(I)
%     dz = z(i) - z(i-1);
%     tmp = tmp + mu(i-1) * dz;
%     I(i) = k*mu(i) .* exp( -2 * tmp ) .* g(i);
%   end

end

