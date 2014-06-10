
function mu = muFitCVX( I, mask, z, z0, zR, eta )

  Md = find( mask, 1, 'last' );
  %Md = Md - 120;

  dz = z(2) - z(1);
  D = makeD_1D(Md, dz);

  h = makeConfocalFunction( z, z0, zR );

  ut = triu( ones(Md-1,Md-1) );
  cvx_begin quiet
    %cvx_precision high

    variable gam(Md)

    minimize ( ...
      0.5 * sum( ...
           ... % Integrate from z to D
           mask(1:Md-1) .* ( ...
            ( 1/2 * I(1:Md-1).* gam(1:Md-1) ./ dz ...
             - ( 1/2 * I(Md) * gam(Md) ./ dz ) ...
             + 1/2 * ut * ( ...
                 I(1:Md-1) ./ dz .* gam(1:Md-1) .* log( h(2:Md)./h(1:Md-1) ) ...
               ) ...
             - ut * I(1:Md-1) ...
             ).^2 ...
           ) ...
      ) ...
      + eta * sum( mask(1:Md) .* abs( D * gam ) ) ...
    )

    subject to
      gam >= 0
      diag(1-mask(1:Md))*gam == 0

  cvx_end

  % Note: we make gamma 0 when mask is 0 because I is 0 at those points

  mu = 1 ./ gam;
  mu( mask == 0 ) = 0;

end
