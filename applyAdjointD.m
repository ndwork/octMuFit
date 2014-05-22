
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


  out = zeros(size(y));
  switch nDimsy
    case 1
      out(1) = -y(1);
      out(end) = y(end-1);
      out(2:end-1) = y(1:end-2) - y(2:end-1);
    case 2
      out(1,:) = -y(1,:);
      out(end,:) = y(end-1,:);
      out(2:end-1,:) = y(1:end-2,:) - y(2:end-1,:);
    case 3
      out(1,:,:) = -y(1,:,:);
      out(end,:,:) = y(end-1,:,:);
      out(2:end-1,:,:) = y(1:end-2,:,:) - y(2:end-1,:,:);
    otherwise
      errror('Improper size of y');
  end
  out = out ./ delta;
  
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






