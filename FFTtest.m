function FFTtest

% Eample use of fftn
% Must first sample the periodic data in a nonrepetative way
% Then do fft of the data to transform

clear
clc

N1 = 10;
N2 = 11;
N3 = 12;
x1vec = linspace(0,1-1/N1,N1);
x2vec = linspace(0,1-1/N2,N2);
x3vec = linspace(0,1-1/N3,N3);

for k = 1:N1
    for l = 1:N2
        for m = 1:N3
            A(k,l,m) = f(x1vec(k),x2vec(l),x3vec(m));
        end
    end
end

fftn(A);

end

% Frequency of 1,1,1
function y = f(x1,x2,x3)

y = exp(2*pi*1i*(x1+x2+x3));

end