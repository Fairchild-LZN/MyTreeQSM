function [V,W] = orthonormal_vectors(U)

% Generate vectors V and W that are unit vectors orthogonal to themselves 
% and to the input vector U

V = rand(3,1);
% 向量叉积
V = cross_product(V,U);
% 当v不是单位向量
while norm(V) == 0
    V = rand(3,1);
    V = cross_product(V,U);
end
W = cross_product(V,U);
% 将w和v变成单位向量
W = W/norm(W);
V = V/norm(V);
if size(V,2) > 1
    V = V';
end
if size(W,2) > 1
    W = W';
end