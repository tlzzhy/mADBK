function [IT,CPU] = RaBK_a_paved(A, repeat)
% ==============================
 % RaBK_a_paved method
 % Input：
 % A:the coefficient matrix
 % repeat: repeat = 10
 % Output:IT：number of iterations, CPU:CPU time
% ==============================
rng(99); 
[m,n]=size(A);
%A = randn(m,n);
temp = sum(A.^2, 2); 
A = full(A);
A = A./sqrt(temp); % normalized
x_true = randn(n,1);
b = A*x_true; 
b_norm = norm(b)^2;
x0 = zeros(n,1);    
sc = zeros(1,repeat);
tc = zeros(1,repeat);
numPartitions = ceil(norm(A, 2)^2);  % Blocksize

% calculate the size of each partition
partitionSize = floor(m / numPartitions);
partitions = cell(numPartitions, 1);
% create partitions
 for i = 1:numPartitions
     startIdx = floor((i - 1) * m / numPartitions) + 1;  % calculate the starting index
     if i < numPartitions
         endIdx = floor(i * m / numPartitions);  % calculate the ending index, excluding the last partition
     else
         endIdx = m;  % the last partition includes all remaining rows
     end
     partitions{i} = startIdx:endIdx;  % store the row index range for each partition
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
        p = randi(numPartitions); % randomly select a partition
        Jk = partitions{p};
        alphak_up = 0;
        dk = zeros(n,1);
        Nk = length(Jk);
        for j = 1:Nk
            i = Jk(j);
            Ai = A(i,:);
            dk = dk + rk(i)* Ai'/temp(i);
            alphak_up = alphak_up + omega * rk(i)^2 / temp(i);
        end
        dk = omega * dk;
        Lk = alphak_up / norm(dk)^2;   
        alphak = 1.95 * Lk; % calculate the step size
       
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