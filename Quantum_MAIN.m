function output = Quantum_MAIN
gui
format long;
clear;
clc;

global L hbar omega m

% Implements variational method for a 1-dimensional problem

L = 2e-10;                                  % length of domain, m
hbar = 1.054571e-34;                        % J*s
omega = 5.63212e14;                         % s^-1
m = 1.62661e-27;                            % Kg

%% Build Potential

t = cputime;

n_bais_vecs = 14;                                     % N + 1 basis functions are used
n_elect = 5;        % number of electrons in the system
n_eng_levels = ceil(n_elect / 2);                                    % number of energy levels retrieved with Lanczos
M = 2 * n_bais_vecs + 1;
x = linspace(-L/2,L/2 - L/M,M);
V = zeros(1,M);

for i = 1:M
    V(i) = harm_pot(x(i));
end

freq = fft(V);                     % Take fft of the potential
freq = circshift(freq,[1,n_bais_vecs]);

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

PW = -n_bais_vecs/2 : n_bais_vecs/2;

% Build kinetic energy matrix
Ham_KE = zeros(n_bais_vecs+1,n_bais_vecs+1);
for k = 1:n_bais_vecs+1
    Ham_KE(k,k) = hbar ^ 2 / 2 / m * L ^ -2 * 4 * pi^2 * PW(k)^2;
end

% Build potential energy matrix
Ham_PE = zeros(n_bais_vecs+1,n_bais_vecs+1);

for i = 1:n_bais_vecs+1
    for j = 1:n_bais_vecs+1
        freqdiff = PW(i) - PW(j);
        index = freqdiff + n_bais_vecs + 1;
        Ham_PE(i,j) = freq(index)/M;
    end
end

Ham = Ham_KE + Ham_PE;
Ham = real(Ham);

FFTtime = cputime - t;

%% Solve SE

t = cputime;

[Vecs, Vals] = eig(Ham);          % returns all eigenvalues and
% eigenvectors
% [Vecs, Vals] = eigs(Ham,n_eng_levels,'SM');       % This does not work...
%Vec = real(Vecs);
%Vals = real(Vals);

%Eigtime = cputime - t

%% Plot Electron Densities
% and compare to analytical solution

eng_lev = 0;                                 % energy level, only ground state wavefunction can be computed anyway
s = 100;                                % sample for graphing
xvec = linspace(-L/2,L/2,s);
xscale = sqrt(hbar/m/omega);            % m

% Analytical Solution
analwf = zeros(1,s);
analpd = zeros(1,s);
for i = 1:s
    analwf(i) = solwf(eng_lev,xvec(i));
    analpd(i) = abs(analwf(i))^2;
end

% Approximate Solution
phimat = zeros(s,n_bais_vecs+1);
for i = 1:s
    for j = 1:n_bais_vecs+1
        phimat(i,j) = phi(PW(j),xvec(i));
    end
end

wfsln = Vecs(:,eng_lev+1);
approxwf = phimat*wfsln;
approxpd = zeros(1,s);
for i = 1:s
    approxpd(i) = abs(approxwf(i))^2;
end

% Non-dimensionalize
xvec = xvec / xscale;
engscale = hbar * omega / 2;
analwf = analwf * sqrt(xscale);
analpd = analpd * xscale;
approxwf = approxwf * sqrt( xscale );
approxpd = approxpd * xscale;

% Flip the wavefunction if it's negative
if approxwf(s/2) < 0
    approxwf = - approxwf;
end

% Compare Energies
E_approx = diag(Vals)';
Eng_lvl = 0:n_bais_vecs;
E_anal = zeros(1,n_bais_vecs+1);

for i = 1:n_bais_vecs+1
    E_anal(i) = Eng(i-1);
end

E_anal = E_anal / engscale;
E_approx = E_approx / engscale;
%GSE_error = E_approx(1) - E_anal(1)

%% Plot results

set(0,'defaultlinelinewidth',1.5)
set(0,'defaultaxeslinewidth',2)

% figure
% plot(PW, real(wfsln), 'o')
% xlabel('Quantum # of PW')
% ylabel('Real Eigenvalue Coeff.')
% ax = gca;
% ax.FontSize = 20;
% 
% figure
% plot(PW, imag(wfsln), 'o')
% xlabel('Quantum # of PW')
% ylabel('Im Eigenvalue Coeff.')
% ax = gca;
% ax.FontSize = 20;
% 
% figure
% plot(xvec, real(analwf), xvec, real(approxwf))
% xlabel('Dimensionles Position')
% ylabel('Real Dimensionless WF')
% legend('Analytical', 'Numerical')
% legend('boxoff')
% ax = gca;
% ax.FontSize = 20;
% 
% figure
% plot(xvec, imag(analwf), xvec, imag(approxwf))
% xlabel('Dimensionles Position')
% ylabel('Im Dimensionless WF')
% legend('Analytical', 'Numerical')
% legend('boxoff')
% ax = gca;
% ax.FontSize = 20;

% figure
% plot(xvec,analpd,xvec,approxpd)
% xlabel('Dimensionles Position')
% ylabel('Dimensionless PDF')
% legend('Analytical', 'Numerical')
% legend('boxoff')
% ax = gca;
% ax.FontSize = 20;

output.a = xvec;
output.b = analpd;
output.c = approxpd;

% figure
% plot(xvec,real(approxwf))
% xlabel('Dimensionles Position')
% ylabel('Real Dimensionless WF')
% ax = gca;
% ax.FontSize = 20;
% 
% figure
% plot(xvec,imag(approxwf))
% xlabel('Dimensionles Position')
% ylabel('IM Dimensionless WF')
% ax = gca;
% ax.FontSize = 20;

% figure
% plot(Eng_lvl,E_anal,Eng_lvl, E_approx)
% xlabel('Energy Level')
% ylabel('Dimensionless Energy')
% legend('Analytical', 'Approximate')

end



%% Harmonic potential

function V = harm_pot(x)

global omega L m


V = 0.5 * m * omega^2 * x ^2;

end

%% Plane wave

function wf = phi(n,x)

global L

xp = x/L + 1/2;

wf = 1/sqrt(L) * exp( 2 * pi * 1i * n * xp);

end

%% Analytical Wavefunction Solution
% n - energy level, starting at 0
function wf = solwf(n,x)

global L hbar omega m

alpha = m * omega / hbar;
y = sqrt(alpha) * x;

wf = (alpha / pi)^0.25 / sqrt(2^n * factorial(n)) * hermite(n,y) * exp(-y^2/2);                            

end

%% Analitical Energy
% n - energy level, starting at 0
function E = Eng(n)

global hbar omega

E = hbar * omega * (n + 0.5);

end



