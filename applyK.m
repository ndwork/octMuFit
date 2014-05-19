
function [yA, yDz, yDx] = applyK( gamma, mask, I, z, dz, dx, z0, zR )

  yA = applyA( gamma, mask, I, dz, z, z0, zR );

  yDz = applyD( gamma, mask, dz, 1 );
  yDx = applyD( gamma, mask, dx, 2 );

end

