function ChemEProject

clc
clear

tic

% Define global variables
global A omega b1 b2 b3 B i1 j1 k1 i2 j2 j3

% Define cell
A = [1 0 0; 0 1 0; 0 0 1]';                  % lengths in m, defined by a1,a2,a3
omega = det(A);                             % Volume of periodic cell

% Compute reciprocal lattice vectors            
b1 = 2 * pi * cross(A(:,2),A(:,3)) / omega;
b2 = 2 * pi * cross(A(:,3),A(:,1)) / omega;
b3 = 2 * pi * cross(A(:,1),A(:,2)) / omega;
B = [b1,b2,b3];

%% Build the plane-wave basis set

E_cut = 170;                                   % Energy cutoff for plane wave basis set

i_max = sqrt(2 * E_cut) / norm(b1,2);
j_max = sqrt(2 * E_cut) / norm(b2,2);
k_max = sqrt(2 * E_cut) / norm(b3,2);

index = 1;
for i = 1:i_max
    for j = 1:j_max
        for k = 1:k_max
            KE_test = 0.5 * (norm(B * [i;j;k],2)) ^ 2;
            if norm(B * [i;j;k],2) < sqrt(2 * E_cut)                         % If the wavefunction with these quantum numbers is below the KE threshold, add it to the basis
                ivec(index) = i;
                jvec(index) = j;
                kvec(index) = k;
                KE(index) = KE_test;                                        % store the kinetic energy value you just calculated
                index = index + 1;
            end
        end
    end
end


PW = [ivec; jvec; kvec]';                           % matrix with the indecies of all the plane waves
index = index - 1                                   % number of plane waves being used, depends only on energy cutoff and lattice size

%% Make Hamiltonian Matrix
Ham_KE = diag(KE);
Ham_PE = zeros(index,index);

for i = 1:index
    for j = 1:index
        Ham_PE(i,j) = Vbraket([PW(i,1), PW(i,2), PW(i,3)], [PW(j,1), PW(j,2), PW(j,3)]);
    end
end
 
% Add kinetic and potential energy contributions to the 
Ham = Ham_KE + Ham_PE

% Find eigenvalues and eigenvectors of A
[V,D] = eigs(Ham);
V                                           % eigenvectors
diag(D)                                     % eigenvalues

toc

end

% Evaluates <phi|V|phi>
function PE = Vbraket(ijk1,ijk2)
global i1 j1 k1 i2 j2 k2
i1 = ijk1(1);
j1 = ijk1(2);
k1 = ijk1(3);
i2 = ijk2(1);
j2 = ijk2(2);
k2 = ijk2(3);
PE = triplequad('V',0,1,0,1,0,1);
end

% Converts fractional to scalar coordinates
function r = frac_to_real(f)
global A
r = A * f;
end