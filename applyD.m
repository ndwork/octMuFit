function out = applyD(gamma, delta, dim)

  nDimsGamma = sum( size(gamma) > 1 );

  if (dim > nDimsGamma)
    error('Dimension argument is too large')
  end


  switch dim
    case 1 %Applying D in the z direction
    case 2 %Applying D in the x direction
      gamma = gamma';
    case 3 %Appyling D in the y direction
      error('Applying D in the third dimension is not supported yet')
    otherwise
      error('Bad dimension argument')
  end

  switch nDimsGamma
    case 1
      out = applyD1D(gamma, delta);
    case 2
      out = applyD2D(gamma, delta);
    case 3
      out = applyD3D(gamma, delta);
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

function out = applyD3D(gamma, delta)
  [M N P] = size(gamma);
  out = zeros(M, N, P);
  for i = 1:N
    for j = 1:M
      out(:, i, j) = applyD1D(gamma(:, i, j), delta);
    end
  end
end

function out = applyD2D(gamma, delta)
  [M N] = size(gamma);
  out = zeros(M, N);
  for i = 1:N
    out(:,i) = applyD1D(gamma(:,i), delta);
  end
end

function out  = applyD1D(gamma, delta)
  out = zeros( numel(gamma), 1 );
  out(1:end-1) = ( gamma(2:end)-gamma(1:end-1) ) ./ delta;
end

