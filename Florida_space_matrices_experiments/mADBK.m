function [IT,CPU] = mADBK(A,beta, repeat)
% ==============================
 % mADBK method
 % Input：
 % A:the coefficient matrix
 % repeat: represents the number of iterations the program is run, which is set to 10 in our paper.
 % beta: In our paper, the momentum parameter beta is set to 0.5.
 % Output:IT：number of iterations, CPU:CPU time
% ==============================
rng(99);
[m,n] = size(A);
%A = randn(m,n);
x_true = randn(n,1); 
b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);
x1 = x0;  
sc = zeros(1,repeat);
tc = zeros(1,repeat);
for r = 1:repeat
    tic;
    p = x0;
    a = x1;
    for k = 1:100000
        rk = b-A*a; 
        error = norm(rk)^2;
        rr = error/b_norm;
        if rr < 1e-6
            tc(r) = toc;
            sc(r) = k; 
            break;
        end
        error_m = error/m;
        Uk = find(rk.^2 >= error_m);
        Nk = length(Uk);
        eta_k = zeros(m,1);
        for ik = 1:Nk
            i = Uk(ik);
            eta_k(i) = rk(i);
        end
        eta = A'*eta_k;
        c= a+norm(eta_k)^2/norm(eta)^2*(eta) + beta * (a - p);
        p = a;
        a = c;
    end
end
CPU = mean(tc);  
IT = mean(sc); 