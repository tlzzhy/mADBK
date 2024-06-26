function [IT,CPU] = RaBK_e_paved(m,n, repeat)
% ==============================
 % RaBK_e_paved method
 % Input：
 % m:the number of rows of the coefficient matrix
 % n:the number of columns in the coefficient matrix
 % repeat: repeat = 10
 % Output:IT：number of iterations, CPU:CPU time
% ==============================
rng(99); 
A = randn(m,n);
temp = sum(A.^2, 2); 
A = A./sqrt(temp);  % normalized
x_true = randn(n,1);
b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);    
sc = zeros(1,repeat);
tc = zeros(1,repeat);
numPartitions = ceil(norm(A, 2)^2);   % Blocksize

% calculate the size of each partition
partitionSize = floor(m / numPartitions);
partitions = cell(numPartitions, 1);
lambdaMaxPerPartition = zeros(numPartitions, 1);  % store the maximum eigenvalue of each partition
% create partitions and calculate the maximum eigenvalue of each partition
 for i = 1:numPartitions
     startIdx = floor((i - 1) * m / numPartitions) + 1;  % calculate the starting index
     if i < numPartitions
         endIdx = floor(i * m / numPartitions);   % calculate the ending index, excluding the last partition
     else
         endIdx = m;   % the last partition includes all remaining rows
     end
     partitions{i} = startIdx:endIdx;   % store the row index range for each partition
     A_sub = A(startIdx:endIdx, :);
     lambdaMaxPerPartition(i) = max(eig(A_sub' * A_sub));
end

omega = 1 / partitionSize;
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
        % randomly select a partition
        p = randi(numPartitions);
        Jk = partitions{p};
        lambda_max = lambdaMaxPerPartition(p);
        alphak = 1.95 * partitionSize / lambda_max;  % calculate the step size based on the current partition
        Nk = length(Jk);
        update = zeros(n, 1);
        for idx = 1:Nk
            ik = Jk(idx);
            ai = A(ik, :);
            update = update + omega*rk(ik) / temp(ik) * ai';
        end
        xn = xk + alphak * update;
    end
end
CPU = mean(tc);  
IT = mean(sc);  

