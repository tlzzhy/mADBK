function [IT,CPU] = GABK(m,n, repeat, zeta)
% ==============================
 % GABK method
 % Input：
 % m:the number of rows of the coefficient matrix
 % n:the number of columns in the coefficient matrix
 % repeat: repeat = 10
 % zeta: The value of zeta is 0.2, which parameters in the GABK algorithm
 % Output:IT：number of iterations, CPU:CPU time
% ==============================
rng(99);
A = randn(m,n);
%[m,n] = size(A);
temp = sum(A.^2, 2); 
x_true = randn(n,1);
b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);    
sc = zeros(1,repeat);
tc = zeros(1,repeat);
for r = 1:repeat
    xk = x0;
    tic;
    for k = 1:100000
        rk = b-A*xk; 
        rr = norm(rk)^2/b_norm;
        if rr < 1e-6
            tc(r) = toc;
            sc(r) = k;  
            break
        end
        gamma_k = rk.^2./temp;
        ek = zeta * max(gamma_k);
        Jk = find(gamma_k >= ek);
        Nk = length(Jk);
        wk = 1/Nk;
        alphak_up = wk * sum(gamma_k(Jk));
        dk = zeros(n,1);
        for j = 1:Nk
            i = Jk(j);
            Ai = A(i,:);
            dk = dk + rk(i)* Ai'/temp(i);
        end
        dk = wk * dk;
        alphak = alphak_up / norm(dk)^2;
        xk = xk + alphak * dk;
    end
end
CPU = mean(tc);  % 求平均迭代时间
IT = mean(sc);  % 求平均迭代步数