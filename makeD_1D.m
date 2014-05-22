
function D = makeD_1D(M, dz)

  % we use symmetric boundary conditions for D
  % this is part of the definition of TV regularization

  D = -eye(M);
  D(M,M) = 0;
  for i=1:M-1
    D(i,i+1) = 1;
  end

  D = D ./ dz;

end
