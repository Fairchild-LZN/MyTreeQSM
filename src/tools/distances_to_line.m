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

% Q内每个点坐标-LinePoint
A = mat_vec_subtraction(Q,LinePoint);
% h内元素和为0
h = A*LineDirec;

% 复制LineDirec
% 复制个数为Q的行数
B = repmat(LineDirec',length(Q(:,1)),1);
B = [h.*B(:,1) h.*B(:,2) h.*B(:,3)];
V = A-B;

d = sqrt(sum(V.*V,2));