function plot_CPU_IT_beta(m,n,repeat,beta)
% ==============================
 % Plot the relationship among the number of iterations, CPU time, and various momentum parameters beta for the mADBK method.
 % Input：
 % m:the number of rows of the coefficient matrix
 % n:the number of columns in the coefficient matrix
 % repeat: The value of repeat is 1
 % beta: In our paper, beta = linspace(0,1,21)
% ==============================
[IT_total_1,CPU_total_1] = beta_plot_mADBK(m,n,repeat,beta);
figure(1)
plot(beta,IT_total_1,'--ob', 'LineWidth',1.0,'MarkerSize',10,'MarkerFaceColor','r');
xlabel('$\beta$','interpreter','latex');
ylabel('IT');
legend('mADBK')
figure(2)
plot(beta,CPU_total_1,'--ob','LineWidth',1.0,'MarkerSize',10,'MarkerFaceColor','r');
xlabel('$\beta$','interpreter','latex');
ylabel('CPU');
legend('mADBK')
%ylabel('$RR$','interpreter','latex');
% ylabel('$\sin(\alpha)$','interpreter','latex', 'FontSize', 18);
end

function [IT_total,CPU_total] = beta_plot_mADBK(m,n,repeat,beta)
rng(99);
%[m,n] = size(A);
A = randn(m,n);
x_true = randn(n,1); b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);
x1 = x0;  
IT_total = zeros(1,length(beta));
CPU_total = zeros(1,length(beta));
for th = 1:length(beta)
    for r = 1:repeat 
    p = x0;
    a = x1;
    tic;
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
        c= a+norm(eta_k)^2/norm(eta)^2*(eta) + beta(th) * (a - p);
        p = a;
        a = c;
    end
    CPU_total(th) = mean(tc);  % 求平均迭代时间
    IT_total(th) = mean(sc);  % 求平均迭代步数
    end
end
end