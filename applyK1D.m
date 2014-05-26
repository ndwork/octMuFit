
function [yA, yDz] = applyK1D( gamma, I, dz, g )

  %M = numel(I);
  %A = makeA_1D( I, dz, g );
  %D = makeD_1D( M, dz );

  %yA = A*gamma;
  %yDz = D*gamma;

  yA = applyA( gamma, I, dz, g );
  yDz = applyD( gamma, dz, 1 );

end

