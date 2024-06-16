function [IT,CPU] = FGBK(m,n, repeat,eta)
% ==============================
 % FGBK method
 % Input：
 % m:the number of rows of the coefficient matrix
 % n:the number of columns in the coefficient matrix
 % repeat:repeat = 10
 % eta: The value of eta is 0.05
 % Output:IT：number of iterations, CPU:CPU time
% ==============================
rng(99); 
%[m,n] = size(A);
A = randn(m,n);
temp = sum(A.^2, 2); 
x_true = randn(n,1);
b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);    
sc = zeros(1,repeat);
tc = zeros(1,repeat);
for r = 1:repeat
    xn = x0;
    tic;
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
        ek = eta * max(rk.^2./temp);
        Uk = find(rk.^2 >= ek * temp);
        Nk = length(Uk);
        eta_k = zeros(m,1);
        for ik = 1:Nk
            i = Uk(ik);
            eta_k(i) = rk(i);
        end
        yk = A'*eta_k;
        xn= xk+norm(eta_k)^2/norm(yk)^2*yk;
    end
end
CPU = mean(tc); 
IT = mean(sc); 