function [Q,H] = LanczosMN( A , k )

for i = 1:k
    
end


% Iteration with j=0 
r0 = rand(D,1); 
b0 = sqrt(r0'*r0); 
q1 = r0/b0;             % random matrix with norm 1 
a1 = q1'*A*q1

%Iteration with j=1
r1 = A*q1 - a1*q1
b1 = sqrt(r1'*r1)
q2 = r1/b1;
a2 = q2'*A*q2

%Iteration with j=2
r2 = A*q2 - a2*q2 - b1*q1;
b2 = sqrt(r2'*r2)
q3 = r2/b2
a3 = q3'*A*q3

% Create Matrix Q
Q = [q1 q2 q3]

%Check orthogonality
EYE = Q'*Q
T = Q'*A*Q

end

