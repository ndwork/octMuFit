
function out = applyAdjointD(y, mask, delta, dim)
  % Apply Adjoint of D
  % Input: y  vector/matrix to apply adjoint of D to
  %        mask   vector/matrix indicating 1 where mu is nonzero
  %        dz     pixel resolution
  %        dim    dimension of y to apply adjoint of D to

  nDimsy = ndims(y);
  if(nDimsy ~= ndims(mask))
      error('y and mask size mismatch');
  end

  if (dim > nDimsy)
      error('Dimension argument is too large')
  end

   switch dim
    case 1 %Applying AdjointD in the z direction
    case 2 %Applying AdjointD in the x direction
      y = y';
      mask = mask';
    case 3 %Appyling AdjointD in the y direction
      error('Applying AdjointD in the third dimension is not supported yet')
    otherwise
      error('Bad dimension argument')
  end


  switch nDimsy
    case 1
      out = applyAdjointD1D(y, mask, delta);
    case 2
      out = applyAdjointD2D(y, mask, delta);
    case 3
      out = applyAdjointD3D(y, mask, delta);
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


function out = applyAdjointD1D(y, mask, delta)

  M = length(y);
  D = makeD_1D(M, mask, delta);
  out = ((D')*y);

end

function out = applyAdjointD2D(y, mask, delta)
  [M N] = size(y);
  out = zeros(M, N);
  
  for i = 1:N
    out(:,i) = applyAdjointD1D(y(:,i), mask(:,i), delta);
  end
end

function out = applyAdjointD3D(y, mask, delta)
    [M N P] = size(y);
    out = zeros(M, N, P);
    for i = 1:N
        for j = 1:M
            out(:, i, j) = applyAdjointD1D(y(:, i, j), mask(:,i, j), delta);
        end
    end
end





