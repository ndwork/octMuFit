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

  out = zeros( size(gamma) );
  switch nDimsGamma
    case 1
      out(1:end-1) = ( gamma(2:end)-gamma(1:end-1) );
    case 2
      out(1:end-1,:) = ( gamma(2:end,:) - gamma(1:end-1,:) );
    case 3
      out(1:end-1,:,:) = ( gamma(2:end,:,:) - gamma(1:end-1,:,:) );
    otherwise
      error('Inproper size of gamma');            
  end
  out = out ./ delta;

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


