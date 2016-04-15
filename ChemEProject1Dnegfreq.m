function ChemEProject1D

format long

clear
clc

global L hbar omega m

% Implements variational method for a 1-dimensional problem

L = 2e-10;                                  % length of domain, m
hbar = 1.054571e-34;                        % J*s
omega = 5.63212e14;                         % s^-1
m = 1.62661e-27;                            % Kg

%% Build Potential

t = cputime;

N = 4;
M = 8 * N + 1;
x = linspace(-L/2,L/2 - L/M,M);
V = zeros(1,M);

for i = 1:M
    V(i) = potential(x(i));
end

freq = fft(V);                     % Take fft of the potential
freq = circshift(freq,[1,4*N]);

%Look at the potential and it's Fourier transform

%figure
%plot(x,V)

% figure
% plot(-N:N,real(freq))
% xlabel('Frequency')
% ylabel('Real part of value')
% 
% figure
% plot(-N:N,imag(freq))
% xlabel('Frequency')
% ylabel('Imaginary part of value')

%% Build Hamiltonian Matrix

PW = -N/2:N/2;
PW = PW([1:N/2, N/2+2:end]);

% Build kinetic energy matrix
Ham_KE = zeros(N,N);
for k = 1:N
    Ham_KE(k,k) = hbar ^ 2 / 2 / m * L ^ -2 * 4 * pi^2 * PW(k)^2;
end

% Build potential energy matrix
Ham_PE = zeros(N,N);
freq
for i = 1:N
    for j = 1:N
        freqdiff = PW(i) - PW(j);
        index = freqdiff + 4*N + 1;
        Ham_PE(i,j) = freq(index)/M;
    end
end

% % Build potential energy matrix
% Ham_PE = zeros(N-1,N-1);
% 
% for i = 1:N-1
%     for j = 1:N-1
%         if j >= i
%             Ham_PE(i,j) = freq(j-i+1)/N;
%         else
%             Ham_PE(i,j) = freq(i-j+1)'/N;
%         end
%     end
% end

Ham = Ham_KE + Ham_PE;

FFTtime = cputime - t

%% Solve SE

t = cputime;

[Vecs, Vals] = eig(Ham);

Eigtime = cputime - t

% Compare energies 
%avg_pot = freq(1)/N
%approx_GSeng = min(diag(Vals))
%GSeng = Eng(0)
%Eng(1)             % look at coefficients from variational method


%% Plot Electron Densities
% and compare to analytical solution

el = 0;                                 % energy level
s = 100;                                % sample for graphing
xvec = linspace(-L/2,L/2,s);
xscale = sqrt(hbar/m/omega)            % m

% Analytical Solution
analwf = zeros(1,s);
analpd = zeros(1,s);
for i = 1:s
    analwf(i) = solwf(el,xvec(i));
    analpd(i) = abs(analwf(i))^2;
end


% xvec = linspace(-L/2,L/2,s);
% yvec = zeros(1,s);
% for i = 1:s
%     yvec(i) = phi(1,xvec(i)) + phi(-1,xvec(i));
% end
% 
% figure
% plot(xvec,yvec)


% Approximate Solution
phimat = zeros(s,N);
for i = 1:s
    for j = 1:N
        phimat(i,j) = phi(PW(j),xvec(i));
    end
end
%Vecs(:,el+1)
wfsln = Vecs(:,el+1);

% %rofl
% for i = N/2+1:N
%     wfsln(i) = -wfsln(i);
% end

approxwf = phimat*wfsln;

% Non-dimensionalize
xvec = xvec / xscale;
analwf = analwf * sqrt(xscale);
approxwf = approxwf * sqrt( xscale );

figure
plot(PW, real(wfsln))
xlabel('Quantum # of PW')
ylabel('real part of Eigenvalue Coefficient')

figure
plot(PW, imag(wfsln))
xlabel('Quantum # of PW')
ylabel('imaginary part ofEigenvalue Coefficient')

figure
plot(xvec,analwf,xvec, approxwf)
xlabel('Dimensionles Position')
ylabel('Dimensionless Wavefunction')
legend('Analytical', 'Approximate')

a = real(approxwf);
b = imag(approxwf);
figure
plot(a,b)
xlabel('Real Part of WF')
ylabel('Imaginary Part of WF')

figure
plot(xvec,a)
xlabel('Position')
ylabel('Real Part of Wavefunction')

figure
plot(xvec,b)
xlabel('Position')
ylabel('Imaginary Part of Wavefunction')

end

function V = potential(x)

global omega L m


V = 0.5 * m * omega^2 * x ^2;

end

function index = freqtoindex(f,M)
    if f >= 0
        index = f + 1;
    else
        index = f + 1 + M;
    end
end

function wf = phi(n,x)

global L

xp = x/L + 1/2;

wf = 1/sqrt(L) * exp( 2 * pi * 1i * n * xp);
wf = real(wf);

end

% Analytical Wavefunction Solution
function wf = solwf(n,x)

global L hbar omega m

alpha = m * omega / hbar;
y = sqrt(alpha) * x;

wf = (alpha / pi)^0.25 / sqrt(2^n * factorial(n)) * hermite(n,y) * exp(-y^2/2);                            

end

% Analitical Energy
function E = Eng(n)

global hbar omega

E = hbar * omega * (n + 0.5);

end



