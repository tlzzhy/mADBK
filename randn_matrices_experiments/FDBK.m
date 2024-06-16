function [IT,CPU] = FDBK(m,n, repeat)
% ==============================
 % FDBK method
 % Input：
 % m:the number of rows of the coefficient matrix
 % n:the number of columns in the coefficient matrix
 % repeat: repeat = 10
 % Output:IT：number of iterations, CPU:CPU time
% ==============================
rng(99);
%[m,n] = size(A);
A = randn(m,n);
temp = sum(A.^2, 2);
Af = sum(temp);
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
        ek = 0.5*(1 / error * max(rk.^2./temp) + 1 /Af); % calculate epsilon_k
        Uk = find(rk.^2 >= ek * error * temp); % calculate the index set 
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