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

N = 10001;
x = linspace(-L/2,L/2 - L/N,N);
V = zeros(1,N);

for i = 1:N
    V(i) = potential(x(i));
end

freq = fft(V);                     % Take fft of the potential


% Look at the potential and it's Fourier transform
% figure
% plot(x,V)
% 
% figure
% plot(0:N-1,freq)

%% Build Hamiltonian Matrix

% Build kinetic energy matrix
Ham_KE = zeros(N-1,N-1);
for k = 1:N-1
    Ham_KE(k,k) = hbar ^ 2 / 2 / m * L ^ -2 * 4 * pi^2 * k^2;
end

% Build potential energy matrix
Ham_PE = zeros(N-1,N-1);

for i = 1:N-1
    for j = 1:N-1
        if j >= i
            Ham_PE(i,j) = freq(j-i+1)/N;
        else
            Ham_PE(i,j) = freq(i-j+1)'/N;
        end
    end
end

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
xscale = sqrt(hbar/m/omega);            % m

% Analytical Solution
analwf = zeros(1,s);
analpd = zeros(1,s);
for i = 1:s
    analwf(i) = solwf(el,xvec(i));
    analpd(i) = abs(analwf(i))^2;
end


% Approximate Solution
phimat = zeros(s,N-1);
for i = 1:s
    for j = 1:N-1
        phimat(i,j) = phi(j,xvec(i));
    end
end
%Vecs(:,el+1)
wfsln = Vecs(:,el+1);
approxwf = phimat*wfsln;

% Non-dimensionalize
xvec = xvec / xscale;
analwf = analwf * sqrt(xscale);
approxwf = approxwf * sqrt( xscale );



figure
plot(xvec,analwf,xvec, approxwf)
xlabel('Dimensionles Position')
ylabel('Dimensionless Wavefunction')
legend('Analytical', 'Approximate')

a = real(approxwf);
b = imag(approxwf);
figure
plot(a,b)

figure
plot(xvec,b)

end

function V = potential(x)

global omega L m

fidgefac1 = 4;
fudgefac2 = 0.4;
fidgefac1 = 1;
fudgefac2 = 1;

V = fudgefac2 * 0.5 * m * omega^2 * (x*fidgefac1) ^2;

end

function wf = phi(n,x)

global L

xp = x/L + 1/2;

wf = 1/sqrt(L) * exp( 2 * pi * 1i * n * xp);

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



