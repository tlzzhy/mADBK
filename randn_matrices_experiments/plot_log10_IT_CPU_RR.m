function plot_log10_IT_CPU_RR(m, n,zeta, eta, beta_mADBK)
% ==============================
 % Plot the relationship between the relative residual RR and both the number of iterations and CPU time for the corresponding method.
 % Input：
 % m:the number of rows of the coefficient matrix
 % n:the number of columns in the coefficient matrix
 % zeta: the value of zeta is 0.2
 % eta: the value of eta is 0.05
 % beta_mADBK: the value of beta_mADBK is 0.5
% ==============================
[IT1, CPU1, FDBK_RR] = FDBK(m,n);
[IT2, CPU2, GABK_RR] = GABK(m,n, zeta);
[IT3, CPU3, FGBK_RR] = FGBK(m,n,eta);
[IT4, CPU4, RSHK_RR] = RSHK(m, n);
[IT5, CPU5, ADBK_RR] = ADBK(m,n);
[IT6, CPU6, mADBK_RR] = mADBK(m,n, beta_mADBK);
figure(1)
semilogy(IT1,FDBK_RR,'--xy',IT2,GABK_RR,'--om', IT3,FGBK_RR,'--<k',IT4,RSHK_RR,'--hr', IT5, ADBK_RR, '-->g',IT6,mADBK_RR,'--*c','LineWidth',1.0,'MarkerSize',8);
legend('FDBK','GABK','FGBK','RSHK', 'ADBK', 'mADBK(0.5)');
grid on;
legend('boxoff');
xlabel('IT');
ylabel('$RR$','interpreter','latex');

figure(2)
semilogy(CPU1,FDBK_RR,'--xy',CPU2,GABK_RR,'--om', CPU3,FGBK_RR,'--<k',CPU4,RSHK_RR,'--hr', CPU5, ADBK_RR, '-->g',CPU6,mADBK_RR,'--*c','LineWidth',0.8,'MarkerSize',6);
legend('FDBK','GABK','FGBK','RSHK', 'ADBK', 'mADBK(0.5)');
grid on;
legend('boxoff');
xlabel('CPU');
ylabel('$RR$','interpreter','latex');
end

function [IT1,CPU1,FDBK_RR] = FDBK(m,n)
rng(99);
A = randn(m,n);
temp = sum(A.^2, 2); 
Af = sum(temp);
x_true = randn(n,1); b = A*x_true; 
xn = zeros(n,1); 
tic;
for k = 1:200000
    xk = xn;
    rk = b-A*xk; 
    error = norm(rk)^2; 
    rr = error/norm(b)^2;
    FDBK_RR(k) = rr;
    IT1(k) = k;
    CPU1(k) = toc;
    if rr < 1e-6 
        break;
    end
    ek = 0.5*(1 / error * max(rk.^2./temp) + 1 /Af);
    Uk = find(rk.^2 >= ek * error * temp);
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

function [IT2,CPU2,GABK_RR] = GABK(m,n, zeta)
rng(99);
A = randn(m,n);
temp = sum(A.^2, 2);
x_true = randn(n,1); b = A*x_true;
xk = zeros(n,1); 
tic
for k = 1:20000
    rk = b-A*xk; 
    rr = norm(rk)^2/norm(b)^2;
    GABK_RR(k) = rr;
    IT2(k) = k;
    CPU2(k) = toc;
    if rr < 1e-6
        break;
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
        dk = dk + rk(i)* Ai'/norm(Ai)^2;
    end
    dk = wk * dk;
    alphak = alphak_up / norm(dk)^2;
    xk = xk + alphak * dk;
end
end

function [IT3, CPU3,FGBK_RR] = FGBK(m,n,eta)
rng(99);
A = randn(m,n);
temp = sum(A.^2, 2);
x_true = randn(n,1);
b = A*x_true; 
xn = zeros(n,1);  
tic;
for k = 1:200000
    xk = xn;
    rk = b-A*xk; 
    error = norm(rk)^2; 
    rr = error/norm(b)^2;
    FGBK_RR(k) = rr;
    IT3(k) = k;
    CPU3(k) = toc;
    if rr < 1e-6
        break;
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

function [IT4,CPU4,RSHK_RR] = RSHK(m, n)
rng(99); 
A = randn(m,n);
x_true = randn(n,1);
b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);    
xn = x0;
tic;
for k = 1:100000
    xk = xn;
    eta_k = b-A*xk;
    eta_k_norm2 = norm(eta_k)^2;
    rr = eta_k_norm2 /b_norm;
    RSHK_RR(k) = rr;
    IT4(k) = k;
    CPU4(k) = toc;
    if rr < 1e-6
        break
    end
    yk = A'*eta_k;
    xn = xk + eta_k_norm2 /norm(yk)^2*yk;
end
end

function [IT5,CPU5,ADBK_RR] = ADBK(m,n)
rng(99); 
A = randn(m,n);
x_true = randn(n,1);
b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);    
xn = x0;
tic;
for k = 1:100000
    xk = xn;
    rk = b-A*xk; 
    error = norm(rk)^2;
    rr = error/b_norm;
    ADBK_RR(k) = rr;
    IT5(k) = k;
    CPU5(k) = toc;
    if rr < 1e-6
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
    xn = xk +norm(eta_k)^2/norm(eta)^2*(eta);
end
end

function [IT6,CPU6,mADBK_RR] = mADBK(m,n, beta_mADBK)
rng(99);
A = randn(m,n);
x_true = randn(n,1); b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);
x1 = x0;
p = x0;
a = x1;
tic;
for k = 1:200000
    rk = b-A*a; 
    error = norm(rk)^2;
    rr = error/b_norm;
    mADBK_RR(k) = rr;
    IT6(k)=k;
    CPU6(k) = toc;
    if rr < 1e-6
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
    c= a+norm(eta_k)^2/norm(eta)^2*(eta) + beta_mADBK * (a - p);
    p = a;
    a = c;
end
end