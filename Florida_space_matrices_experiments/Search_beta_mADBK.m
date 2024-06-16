function SCPU = Search_beta_mADBK(A,beta)
% calculate the CPU time for searching the optimal 𝛽
% beta:beta=linspace(0,1,21)
rng(99);
[m,n] = size(A);
%A = randn(m,n);
x_true = randn(n,1); b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);
x1 = x0;  
CPU_total = zeros(1,length(beta));
IT_total = zeros(1,length(beta));
for th = 1:length(beta)
    p = x0;
    a = x1;
    tic;
    for k = 1:100000
        rk = b-A*a; 
        error = norm(rk)^2;
        rr = error/b_norm;
        if rr < 1e-6
            CPU_total(th) = toc;
            IT_total(th) = k;
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
        c= a+norm(eta_k)^2/norm(eta)^2*(eta) + beta(th) * (a - p);
        p = a;
        a = c;
    end
end
disp(CPU_total)
disp(IT_total)
SCPU = sum(CPU_total);
end