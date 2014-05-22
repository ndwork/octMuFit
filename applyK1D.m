
function [yA, yDz] = applyK1D( gamma, I, z, dz, z0, zR )

  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );

  %M = numel(I);
  %A = makeA_1D( I, dz, g );
  %D = makeD_1D( M, dz );

  %yA = A*gamma;
  %yDz = D*gamma;

  yA = applyA( gamma, I, dz, z, z0, zR );
  yDz = applyD( gamma, dz, 1 );

end

