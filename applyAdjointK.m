
function adjKy = applyAdjointK( yA, yDz, yDx, I, ...
  mask, z, dz, dx, z0, zR )

  adjAY = applyAdjointA( yA, mask, I, z, dz, z0, zR );

  adjDzY = applyAdjointD( yDz, mask, dz, 1 );

  adjDxY = applyAdjointD( yDx, mask, dx, 2 );

  adjKy = adjAY + adjDzY + adjDxY;
end
