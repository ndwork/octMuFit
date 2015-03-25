
function muFit = muFit2D_TV( I, z, z0, zR )

  theta = 1.0;    % regularization parameter

  muFit_mVer = muFit2D_mVer( I, z, z0, zR );

  muFit = ROFdenoise(muFit_mVer, theta);

end
