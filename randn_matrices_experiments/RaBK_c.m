function [IT2,CPU2] = RaBK_c(m, n, tau, repeat)
% ==============================
 % RaBK_c method
 % Input：
 % m:the number of rows of the coefficient matrix
 % n:the number of columns in the coefficient matrix
 % tau:the value of tau is 10
 % repeat: repeat = 10
 % Output:IT：number of iterations, CPU:CPU time
% ==============================
rng(99); 
A = randn(m,n);
temp = sum(A.^2, 2); 
x_true = randn(n,1);
b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);    
sc = zeros(1,repeat);
tc = zeros(1,repeat);
omega = 1 / tau;
alpha = 1.95; % step size
for r = 1:repeat
    tic;
    xn = x0;
    for k = 1:100000
        xk = xn;
        rk = b-A*xk; 
        error = norm(rk)^2; 
        rr = error/b_norm;
        if rr < 1e-6
            tc(r) = toc;
            sc(r) = k; 
            break
        end
        Jk = randperm(m, tau);  % select 𝜏 row indices using uniform sampling
        update = zeros(n, 1);
        for idx = 1:tau
            i = Jk(idx);
            ai = A(i, :);
            update = update + omega * rk(i) / temp(i) * ai'; 
        end
        xn = xk + alpha * update;
    end
end
CPU2 = mean(tc); 
IT2 = mean(sc);  

