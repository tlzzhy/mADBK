function [IT,CPU] = ADBK(m,n, repeat)
% ==============================
 % ADBK method
 % Input：
 % m:the number of rows of the coefficient matrix
 % n:the number of columns in the coefficient matrix
 % repeat: repeat = 10
 % Output:IT：number of iterations, CPU:CPU time
% ==============================
rng(99); % set the random seed
%[m,n] = size(A);
A = randn(m,n);
x_true = randn(n,1);
b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1); % initial iterate vector
sc = zeros(1,repeat); % 
tc = zeros(1,repeat);
for r = 1:repeat
    tic;
    xn = x0;
    for k = 1:100000
        xk = xn;
        rk = b-A*xk; 
        error = norm(rk)^2;
        rr = error/b_norm;
        if rr < 1e-6
            tc(r) = toc; % store the CPU time
            sc(r) = k;  % store the number of iterations
            break;
        end
        error_m = error/m;
        Uk = find(rk.^2 >= error_m); % the set of positive integers in step three of the algorithm
        Nk = length(Uk);
        eta_k = zeros(m,1);
        for ik = 1:Nk
            i = Uk(ik);
            eta_k(i) = rk(i);
        end
        eta = A'*eta_k;
        xn = xk +norm(eta_k)^2/norm(eta)^2*(eta);
    end
end
CPU = mean(tc);  % calculate the average CPU time
IT = mean(sc); 