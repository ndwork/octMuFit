
function adjKy = applyAdjointK1D( yA, yDz, I, dz, g )

  %M = numel(I);
  %A = makeA_1D( I, dz, g );
  %D = makeD_1D( M, mask, dz );

  %adjAY = A'*yA;
  %adjDzY = D'*yDz;


  adjAY = applyAdjointA( yA, I, dz, g );
  adjDzY = applyAdjointD( yDz, dz, 1 );

  adjKy = adjAY + adjDzY;
end
