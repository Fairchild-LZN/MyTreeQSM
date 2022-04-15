function [Vector,IndElements] = verticalcat(CellArray)

% CellArray 输入的分支每层layer的集合索引
    
% Vertical concatenation of the given cell-array into a vector.

% 求每层layer包括的集合个数
CellSize = cellfun('length',CellArray); % determine the size of each cell
% 有多少层
nc = max(size(CellArray)); % number of cells
% 初始化一个ncx2的矩阵
IndElements = ones(nc,2); % indexes for elements in each cell
% 计算每层的累计值，保存在第二列
IndElements(:,2) = cumsum(CellSize);
% 将第二列下移一层，copy到第一列中
IndElements(2:end,1) = IndElements(2:end,1)+IndElements(1:end-1,2);
Vector = zeros(sum(CellSize),1); % concatenation of the cell-array into a vector
for j = 1:nc
    Vector(IndElements(j,1):IndElements(j,2)) = CellArray{j};
end