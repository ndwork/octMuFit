
function out = applyAdjointD(y, delta, dim)
  % Apply Adjoint of D
  % Input: y  vector/matrix to apply adjoint of D to
  %        dz     pixel resolution
  %        dim    dimension of y to apply adjoint of D to

  nDimsy = ndims(y);

  if (dim > nDimsy)
      error('Dimension argument is too large')
  end

   switch dim
    case 1 %Applying AdjointD in the z direction
    case 2 %Applying AdjointD in the x direction
      y = y';
    case 3 %Appyling AdjointD in the y direction
      error('Applying AdjointD in the third dimension is not supported yet')
    otherwise
      error('Bad dimension argument')
  end


  switch nDimsy
    case 1
      out = applyAdjointD1D(y, delta);
    case 2
      out = applyAdjointD2D(y, delta);
    case 3
      out = applyAdjointD3D(y, delta);
    otherwise
      errror('Improper size of y');
  end
  
  switch dim
    case 1 %Applying AdjointD in the z direction
    case 2 %Applying AdjointD in the x direction
      out = out';
    case 3 %Appyling AdjointD in the y direction
      error('Applying AdjointD in the third dimension is not supported yet')
    otherwise
      error('Bad dimension argument')
   end

end


function out = applyAdjointD1D(y, delta)

  M = length(y);
  D = makeD_1D(M, delta);
  out = ((D')*y);

end

function out = applyAdjointD2D(y, delta)
  [M N] = size(y);
  out = zeros(M, N);
  
  for i = 1:N
    out(:,i) = applyAdjointD1D(y(:,i), delta);
  end
end

function out = applyAdjointD3D(y, delta)
    [M N P] = size(y);
    out = zeros(M, N, P);
    for i = 1:N
        for j = 1:M
            out(:, i, j) = applyAdjointD1D(y(:, i, j), delta);
        end
    end
end





