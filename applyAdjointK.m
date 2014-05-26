
function adjKy = applyAdjointK( yA, yDz, yDx, I, dz, dx, g )

  adjAY = applyAdjointA(yA, I, dz, g);

  adjDzY = applyAdjointD(yDz, dz, 1);

  adjDxY = applyAdjointD(yDx, dx, 2);

  adjKy = adjAY + adjDzY + adjDxY;
end
