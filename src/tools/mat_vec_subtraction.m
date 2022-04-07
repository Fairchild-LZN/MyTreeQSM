function A = mat_vec_subtraction(A,v)

% Subtracts from each row of the matrix A the vector v.
% If A is (n x m)-matrix, then v needs to be m-vector.

% A-V，每次计算一列
% 这种方法减只需要三次循环，xyz就可以结束
% 如果按行减，需要根据矩阵内的元素个数决定
for i = 1:size(A,2)
    A(:,i) = A(:,i)-v(i);
end