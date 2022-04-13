function [d,V,h,B] = distances_to_line(Q,LineDirec,LinePoint)

% LineDirec方向向量
% LinePoint是底部Base的中心点的中心值

% Calculates the distances of the points, given in the rows of the
% matrix Q, to the line defined by one of its point and its direction.

if size(LineDirec,1) == 1
    LineDirec = LineDirec';
end
% 再次确认已经是 单位向量
LineDirec = LineDirec/norm(LineDirec);

% Q内每个点坐标-LinePoint（Base点）
% Q是base各个集合的中心点
% A就是树的xoy一圈，各个方向朝着中心点的各个向量
A = mat_vec_subtraction(Q,LinePoint);
% h内元素和为0
% LineDirec（Base的法向量），法向量和各个垂直量相乘
h = A*LineDirec;

% 复制LineDirec
% 复制个数为Q的行数，根据base中的集合个数复制
B = repmat(LineDirec',length(Q(:,1)),1);
% 需要了解一下向量相乘的几何含义
B = [h.*B(:,1) h.*B(:,2) h.*B(:,3)];
V = A-B;

d = sqrt(sum(V.*V,2));