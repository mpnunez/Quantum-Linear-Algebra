function ChemEProject2

clc
clear

tic

% Define global variables
%global A omega b1 b2 b3 B

% Define cell
A = [1 0 0; 0 1 0; 0 0 1]';                  % lengths in m, defined by a1,a2,a3
omega = det(A);                             % Volume of periodic cell

% Compute reciprocal lattice vectors            
b1 = 2 * pi * cross(A(:,2),A(:,3)) / omega;
b2 = 2 * pi * cross(A(:,3),A(:,1)) / omega;
b3 = 2 * pi * cross(A(:,1),A(:,2)) / omega;
B = [b1,b2,b3];

%% FFT of the potential
N = 3;
N1 = N; N2 = N; N3 = N;                     % Make it possible to change sampling along different vectors

% Coordinates for sampling
x1vec = linspace(0,1-1/N1,N1);
x2vec = linspace(0,1-1/N2,N2);
x3vec = linspace(0,1-1/N3,N3);

sam = zeros(N1,N2,N3);                      % Multidimensional array with sample points of the potential
PWindices = zeros(N1*N2*N3,4);              

index = 1;
for k = 1:N1
    for l = 1:N2
        for m = 1:N3
            % Sample the point for FFT
            sam(k,l,m) = pot(x1vec(k),x2vec(l),x3vec(m));
            
            % Evaluate PW KE and store quantum numbers
            KE = 0.5 * (norm(B * [k-1;l-1;m-1],2)) ^ 2;
            PWindices(index,:) = [k-1,l-1,m-1,KE];
            index = index + 1;
        end
    end
end

PWindices
fftn(sam);

%% Build the plane-wave basis set



%% Make Hamiltonian Matrix
Ham_KE = diag(PWindices(:,4));

toc

end

% Potential, function of fractional coordinates
% Frequency of 1,2,3
function y = pot(x1,x2,x3)

y = exp(2*pi*1i*(x1+2*x2+3*x3));

end