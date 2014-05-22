
function adjKy = applyAdjointK1D( yA, yDz, I, z, dz, z0, zR )


  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );

  %M = numel(I);
  %A = makeA_1D( I, dz, g );
  %D = makeD_1D( M, mask, dz );

  %adjAY = A'*yA;
  %adjDzY = D'*yDz;


  adjAY = applyAdjointA( yA, I, z, dz, z0, zR );
  adjDzY = applyAdjointD( yDz, dz, 1 );

  adjKy = adjAY + adjDzY;
end
