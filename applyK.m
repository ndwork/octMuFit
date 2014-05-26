
function [yA, yDz, yDx] = applyK( gamma, I, dz, dx, g )

  yA = applyA( gamma, I, dz, g );

  yDz = applyD( gamma, dz, 1 );
  yDx = applyD( gamma, dx, 2 );

end

