
function [yA, yDz, yDx] = applyK( gamma, I, z, dz, dx, z0, zR )

  yA = applyA( gamma, I, dz, z, z0, zR );

  yDz = applyD( gamma, dz, 1 );
  yDx = applyD( gamma, dx, 2 );

end

