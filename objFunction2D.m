function out = objFunction2D(gam,I,mask,dz,dx,z,z0,zR,eta)
  
  b = makeb(I);
  Ag = applyA(gam, I, dz, z, z0, zR);
  Dgz = applyD(gam, dz, 1);
  Dgx = applyD(gam, dx, 2);
  
  out = 0.5 * sum ( mask.*(Ag - b).^2 ) + eta*sum(abs(mask.*Dgz)) + eta*sum(abs(mask.*Dgx));

end
