function [Vector,IndElements] = verticalcat(CellArray)

% Vector返回当前分支按顺序排列的集合索引
% IndElements返回当前分支每层的个数

% Vertical concatenation of the given cell-array into a vector.

% 当前分支每层包括多少个集合
CellSize = cellfun('length',CellArray); % determine the size of each cell
% 当前分支有多少层
nc = max(size(CellArray)); % number of cells
IndElements = ones(nc,2); % indexes for elements in each cell
% 第二列是每层集合个数的累计
IndElements(:,2) = cumsum(CellSize);
% 将第二列向左、向下移动一格，并+1
% 这两步的目的是为了方便后续的赋值
IndElements(2:end,1) = IndElements(2:end,1)+IndElements(1:end-1,2);
% 当前分支一共包括多少个集合
Vector = zeros(sum(CellSize),1); % concatenation of the cell-array into a vector
for j = 1:nc
    % 将相对应的集合索引，赋值给vector
    Vector(IndElements(j,1):IndElements(j,2)) = CellArray{j};
end