function [d,V,h,B] = distances_to_line(Q,LineDirec,LinePoint)

% LineDirec方向向量
% LinePoint是底部Base的中心点的中心值

% Calculates the distances of the points, given in the rows of the
% matrix Q, to the line defined by one of its point and its direction.

% 

if size(LineDirec,1) == 1
    LineDirec = LineDirec';
end
% 再次确认已经是 单位向量
LineDirec = LineDirec/norm(LineDirec);

% Q内每个点坐标-LinePoint（Base点）
% Q是base各个集合的中心点
% A就是树的xoy一圈，中心点朝着，每个点的各个向量
A = mat_vec_subtraction(Q,LinePoint);

% 0417更新，下面这两行可能有问题
% 可能只适用于，A和LineDirec互相垂直的情况
% h内元素和为0
% LineDirec（Base的法向量），法向量和各个垂直量相乘

% 0417更新
% 本质向量相乘，点积，A*B = A * B * cos角度
h = A*LineDirec;

% 复制LineDirec
% 复制个数为Q的行数，根据base中的集合个数复制
B = repmat(LineDirec',length(Q(:,1)),1);
% 需要了解一下向量相乘的几何含义

% 0417更新
% 本质向量相乘，点积，A*B = A * B * cos角度

B = [h.*B(:,1) h.*B(:,2) h.*B(:,3)];

% 0417更新
% 究极本质 V = (A-Center)(1-a^2)
% A是每个点的三维坐标，a是单位方向向量, center是最底部的中心点
V = A-B;

d = sqrt(sum(V.*V,2));

% 0417我觉得我懂了
% 但应该还是不懂