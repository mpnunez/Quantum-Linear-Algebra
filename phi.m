% Value of plane wave, uses fractional coordinates
function phi = phi(ijk,f1,f2,f3)

m = length(f1);
phi = zeros(m,1);

for index = 1:m
    phi(index) = exp(2*pi*complex(0,1) * dot(ijk,[f1(index),f2,f3]));
end

end
