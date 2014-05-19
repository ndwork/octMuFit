function out = applyD(gamma, mask, delta, dim)

  nDimsGamma = sum( size(gamma) > 1 );

  ndimsMask = sum( size(mask) > 1 );
  if(nDimsGamma ~= ndimsMask)
    error('gamma and mask size mismatch');
  end

  if (dim > nDimsGamma)
    error('Dimension argument is too large')
  end


  switch dim
    case 1 %Applying D in the z direction
    case 2 %Applying D in the x direction
      gamma = gamma';
      mask = mask';
    case 3 %Appyling D in the y direction
      error('Applying D in the third dimension is not supported yet')
    otherwise
      error('Bad dimension argument')
  end

  switch nDimsGamma
    case 1
      out = applyD1D(gamma, mask, delta);
    case 2
      out = applyD2D(gamma, mask, delta);
    case 3
      out = applyD3D(gamma, mask, delta);
    otherwise
      error('Inproper size of gamma');            
  end

  switch dim
    case 1 %Applying D in the z direction
    case 2 %Applying D in the x direction
      out = out';
    case 3 %Appyling D in the y direction
      error('Applying D in the third dimension is not supported yet')
    otherwise
      error('Bad dimension argument')
  end

end

function out = applyD3D(gamma, mask, delta)
  [M N P] = size(gamma);
  out = zeros(M, N, P);
  for i = 1:N
    for j = 1:M
      out(:, i, j) = applyD1D(gamma(:, i, j), mask(:,i, j), delta);
    end
  end
end

function out = applyD2D(gamma, mask, delta)
  [M N] = size(gamma);
  out = zeros(M, N);
  for i = 1:N
    out(:,i) = applyD1D(gamma(:,i), mask(:,i), delta);
  end
end

function out  = applyD1D(gamma, mask, delta)
  out = zeros( numel(gamma), 1 );
  out(1:end-1) = ( gamma(2:end)-gamma(1:end-1) ) ./ delta;
  shifted = circshift(mask, [-1, 0]);
	shifted(end,:) = 1;
  tmp = mask & shifted;
  out = out .* tmp;
end

