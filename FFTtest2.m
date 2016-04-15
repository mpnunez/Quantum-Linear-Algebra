function FFTtest2

clear
clc

N = 20;
x = linspace(0,1-1/N,N);
y = zeros(1,N);
for i = 1:N
    if x(i) > 0.4 && x(i) < 0.6
        y(i) = 1;
    end 
end

f = fft(y)
z = ifft(f)
plot(x,y)
%plot(x,z)

end

