function [IT,CPU] = RSHK(A, repeat)
% ==============================
 % RSHK method
 % Input：
 % A:the coefficient matrix
 % repeat: represents the number of iterations the program is run, which is set to 10 in our paper.
 % Output:IT：number of iterations, CPU:CPU time
% ==============================
rng(99); 
[m,n] = size(A);
%A = randn(m,n);
x_true = randn(n,1);
b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);    
% 循环之前赋值
sc = zeros(1,repeat);
tc = zeros(1,repeat);
tic;
for r = 1:repeat
    tic;
    xn = x0;
    for k = 1:100000
        xk = xn;
        eta_k = b-A*xk;
        eta_k_norm2 = norm(eta_k)^2;
        rr = eta_k_norm2 /b_norm;
        if rr < 1e-6
            tc(r) = toc;
            sc(r) = k; 
            break
        end
        yk = A'*eta_k;
        xn = xk + eta_k_norm2 /norm(yk)^2*yk;
    end
end
CPU = mean(tc); 
IT = mean(sc); 
