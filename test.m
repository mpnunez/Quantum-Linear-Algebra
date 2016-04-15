function test

% figure out how to get eigenvector from this

clear
clc

n = 1000;

A = rand(n,n);                 
%A = A + A';
[V,D] = eigs(A);           % matlab procedure
diag(D)


% % Compare exact and approximate eigenvalues and eigenvectors
% eigval          % extract highest magnitude eigenvalues
% eigval_approx = eig(b)
% eigvec
% eigvec_approx = V2(:,1:k)*a                            % approximate eigenvectors of the original matrix


% [sortedValues,sortIndex] = sort(ds(:),'descend');  %# Sort the values in descending order
% maxIndex = sortIndex(1:5)  %# Get a linear index into A of the 5 largest values

%% Test Lanczos

% % Create random Hermetian matrix
% n = 5;
% A = zeros(n,n);
% for i = 1:n
%     A(i,i) = rand;
% end
% for j = 2:n
%     for i = 1:j-1
%         A(i,j) = rand + complex(0,1)*rand;
%         A(j,i) = A(i,j)';
%     end
% end

%% Test Arnoldi code

% n = 20;
% A = rand(n,n);
% b0 = rand(n,1);
% b0 = b0 / norm(b0,2)
% [Q,H] = arnoldi3(A,b0,3)
% 
% size(Q)
% size(H)

end

