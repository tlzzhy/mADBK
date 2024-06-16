function [IT,CPU] = RaBK_a(m, n, tau, repeat)
% ==============================
 % RaBK_a method
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
        Jk = randperm(m, tau); % select 𝜏 row indices using uniform sampling
        gamma_k = rk.^2./temp;
        wk = 1/tau;
        alphak_up = wk * sum(gamma_k(Jk));
        dk = zeros(n,1);
        Nk = length(Jk);
        for j = 1:Nk
            i = Jk(j);
            Ai = A(i,:);
            dk = dk + rk(i)* Ai'/temp(i);
        end
        dk = wk * dk;
        Lk = alphak_up / norm(dk)^2;   
        alphak = 1.95 * Lk;  % calculate the step size
        update = zeros(n, 1);
        for idx = 1:tau
            i = Jk(idx);
            ai = A(i, :);
            update = update + omega * rk(i) / temp(i) * ai'; 
        end
        % 更新解向量 x
        xn = xk + alphak * update;
    end
end
CPU = mean(tc);  % 求平均迭代时间
IT = mean(sc);  % 求平均迭代步数
