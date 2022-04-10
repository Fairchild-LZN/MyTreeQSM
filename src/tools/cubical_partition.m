% This file is part of TREEQSM.
% 
% TREEQSM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TREEQSM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TREEQSM.  If not, see <http://www.gnu.org/licenses/>.

function [Partition,CubeCoord,Info,Cubes] = cubical_partition(P,EL,NE)
% P为输入矩阵
% EL为分割的边长长度 = 0.095 = 0.08 + 0.015     create_inputs.m
% NE为长宽高分割时，两头多出的空白正方形


% Partitions the input point cloud into cubes.
%
% Inputs:
% P             Point cloud, (n_points x 3)-matrix
% EL            Length of the cube edges
% NE            Number of empty edge layers
%
% Outputs:              
% Partition     Point cloud partitioned into cubical cells,
%                   (nx x ny x nz)-cell, where nx,ny,nz are the number
%                   of cubes in x,y,z-directions, respectively
% CC            (n_points x 3)-matrix whose rows are the cube coordinates 
%                   of each point: x,y,z-coordinates
% Info          The minimum coordinate values and number of cubes in each
%                   coordinate direction

if nargin == 2
    NE = 3;
end

% The vertices of the big cube containing P
% 返回P内每一列最小的元素
Min = double(min(P));
% 返回P内每一列最大的元素
Max = double(max(P));

% Number of cubes with edge length "EdgeLength" in the sides 
% of the big cube
% 分割成边长为EL的小正方形，N保存长宽高的正方形个数，总数为三项相乘
% +1 不太理解为什么
% 0327更新，+1的目的是为了让所有都被分割，因为可能会出现不是整数的情况，而多余的部分舍掉也没事的
N = double(ceil((Max-Min)/EL)+2*NE+1);
% 如果正方形总数太多，则扩大分割的边长长度
while 8*N(1)*N(2)*N(3) > 4e9
    EL = 1.1*EL;
    N = double(ceil((Max-Min)/EL)+2*NE+1);
end
% 从左往右，依次为三维最小值、三维分块数、边长、两头多出的正方形数
Info = [Min N EL NE];

% Calculates the cube-coordinates of the points
% 将矩阵P归一化，长宽高的最小值变为0，然后确定每个点，分别位于长宽高中哪儿块正方形内
CubeCoord = floor([P(:,1)-Min(1) P(:,2)-Min(2) P(:,3)-Min(3)]/EL)+NE+1;

% Sorts the points according a lexicographical order
% nx3矩阵 乘 3x1矩阵 = nx1矩阵
% 不理解为什么相乘
LexOrd = [CubeCoord(:,1) CubeCoord(:,2)-1 CubeCoord(:,3)-1]*[1 N(1) N(1)*N(2)]';
CubeCoord = uint16(CubeCoord);
% 给LexOrd升序排序,SortOrd为排序前的位置索引值
[LexOrd,SortOrd] = sort(LexOrd);
SortOrd = uint32(SortOrd);
LexOrd = uint32(LexOrd);

% nargout内置函数,判断函数输出个数
if nargout <= 3
    % Define "Partition"
    % 创建[N(1)xN(2)xN(3)]的矩阵
    Partition = cell(N(1),N(2),N(3));
    np = size(P,1);     % number of points
    p = 1;              % The index of the point under comparison

    % 根据LexOrd计算并排序后的结果
    % 寻找到计算与计算结果相同的点--我称为x
    % 并把这些值的原始位置y--在P矩阵内位于第几个点
    % 找到x中某一个点,位于长宽高正方形的位置z
    % 把y值保存在Partition的z位置处
    while p <= np
        t = 1;
        while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
            t = t+1;
        end
        q = SortOrd(p);
        Partition{CubeCoord(q,1),CubeCoord(q,2),CubeCoord(q,3)} = SortOrd(p:p+t-1);
        p = p+t;
    end

% else代码没看
else
    % Lexord中不同值的个数
    nc = size(unique(LexOrd),1);
    
    % Define "Partition"
    % N(1)xN(2)xN(3)维的全0矩阵，正方形的长宽高三维个数
    Cubes = zeros(N(1),N(2),N(3),'uint32');
    % ncx1维度的元胞数组
    Partition = cell(nc,1);
    np = size(P,1);     % number of points
    p = 1;              % The index of the point under comparison
    c = 0;
    while p <= np
        t = 1;
        % 找到值一样的点坐标
        while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
            t = t+1;
        end
        q = SortOrd(p);
        c = c+1;
        % 将位于P内的索引值，存入Partition中
        Partition{c,1} = SortOrd(p:p+t-1);
        % Cubes中有很多很多的空白，一共只有nc个位置有值
        Cubes(CubeCoord(q,1),CubeCoord(q,2),CubeCoord(q,3)) = c;
        p = p+t;
    end
end