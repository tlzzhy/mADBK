function [IT2,CPU2] = RBK(m,n,block_num, repeat)
% ==============================
 % RBK method
 % Input：
 % m:the number of rows of the coefficient matrix
 % n:the number of columns in the coefficient matrix
 % block_num: the number of blocks into which the coefficient matrix is divided
 % repeat: represents the number of iterations the program is run, which is set to 10 in our paper.
 % Output:IT：number of iterations, CPU:CPU time
% ==============================
rng(99); 
blockSize = m/block_num; % In all our numerical experiments, the blocksize is set to 10.
A = randn(m,n);
[m,n] = size(A);
x_true = randn(n,1);
b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);    
sc = zeros(1,repeat);
tc = zeros(1,repeat);
for r = 1:repeat
    xn = x0;
    tic;
    T = partitionMatrix(m, block_num, blockSize);
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
        tau = T{randi(length(T))};  % randomly select a block
        A_tau = A(tau, :);
        b_tau = b(tau); 
        xn = xk + lsqr(A_tau, b_tau - A_tau * xk, 1e-6, n);
    end
end
CPU2 = mean(tc);  
IT2 = mean(sc); 
end

function T = partitionMatrix(m, block_num, blockSize)
    % initialize the partition array
    T = cell(1, block_num);
    for i = 1:block_num
        startIdx = (i-1) * blockSize + 1;
        endIdx = min(i * blockSize, m);  % ensure the last block does not exceed the index range
        T{i} = startIdx:endIdx;
    end
end

