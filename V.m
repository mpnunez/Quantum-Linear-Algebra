% Potential function used in the calculation
function V = V(f1,f2,f3)

global i1 j1 k1 i2 j2 k2

a = phi([i1, j1, k1],f1,f2,f3);
b = 0*exp( -(f1-0.5).^2 - (f2-0.5).^2 - (f3-0.5).^2 )';
c = phi([i2, j2, k2],f1,f2,f3);

V = a .* b .* c;

end

