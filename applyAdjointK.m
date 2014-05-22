
function adjKy = applyAdjointK( yA, yDz, yDx, I, z, dz, dx, z0, zR )

  adjAY = applyAdjointA(yA, I, z, dz, z0, zR);

  adjDzY = applyAdjointD(yDz, dz, 1);

  adjDxY = applyAdjointD(yDx, dx, 2);

  adjKy = adjAY + adjDzY + adjDxY;
end
