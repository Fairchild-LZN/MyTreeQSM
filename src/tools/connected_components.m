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

function [Components,CompSize] = connected_components(Nei,Sub,MinSize,Fal)

% ---------------------------------------------------------------------
% CONNECTED_COMPONENTS.M      Determines the connected components of cover
%                                   sets using their neighbour-relation
%
% Version 1.1
% Latest update     16 Aug 2017
%
% Copyright (C) 2013-2017 Pasi Raumonen
% ---------------------------------------------------------------------

% Determines connected components of the subset of cover sets defined
% by "Sub" such that each component has at least "MinSize"
% number of cover sets.
%
% Inputs:
% Nei       Neighboring cover sets of each cover set, (n_sets x 1)-cell
% Sub       Subset whose components are determined,
%               length(Sub) < 2 means no subset and thus the whole point cloud
%               "Sub" may be also a vector of cover set indexes in the subset
%               or a logical (n_sets)-vector, where n_sets is the number of
%               all cover sets
% MinSize   Minimum number of cover sets in an acceptable component
% Fal       Logical false vector for the cover sets
%
% Outputs:
% Components    Connected components, (n_comp x 1)-cell
% CompSize      Number of sets in the components, (n_comp x 1)-vector

% Sub是cover_sets生成的集合数
% Sub内的变量若为0,则代表已连接
% Sub内的变量若为1,则代表该集合没有连接

% length(Sub)一般一定大于3个，所以不会进入第一个 if

if length(Sub) <= 3 && ~islogical(Sub) && Sub(1) > 0
    % Very small subset, i.e. at most 3 cover sets
    n = length(Sub);
    if n == 1
        Components = cell(1,1);
        Components{1} = uint32(Sub);
        CompSize = 1;
    elseif n == 2
        I = Nei{Sub(1)} == Sub(2);
        if any(I)
            Components = cell(1,1);
            Components{1} = uint32((Sub));
            CompSize = 1;
        else
            Components = cell(2,1);
            Components{1} = uint32(Sub(1));
            Components{2} = uint32(Sub(2));
            CompSize = [1 1];
        end
    elseif n == 3
        I = Nei{Sub(1)} == Sub(2);
        J = Nei{Sub(1)} == Sub(3);
        K = Nei{Sub(2)} == Sub(3);
        if any(I)+any(J)+any(K) >= 2
            Components = cell(1,1);
            Components{1} = uint32(Sub);
            CompSize = 1;
        elseif any(I)
            Components = cell(2,1);
            Components{1} = uint32(Sub(1:2));
            Components{2} = uint32(Sub(3));
            CompSize = [2 1];
        elseif any(J)
            Components = cell(2,1);
            Components{1} = uint32(Sub([1 3]));
            Components{2} = uint32(Sub(2));
            CompSize = [2 1];
        elseif any(K)
            Components = cell(2,1);
            Components{1} = uint32(Sub(2:3));
            Components{2} = uint32(Sub(1));
            CompSize = [2 1];
        else
            Components = cell(3,1);
            Components{1} = uint32(Sub(1));
            Components{2} = uint32(Sub(2));
            Components{3} = uint32(Sub(3));
            CompSize = [1 1 1];
        end
    end

% 当Sub内全为0时，返回0，有一个不是0，则返回1
elseif any(Sub) || (length(Sub) == 1 && Sub(1) == 0)
    nb = size(Nei,1);
    if nargin == 3
        % 若输入参数为3个,则设Fal全为0;  实际情况下,输入的Fal参数也全为0
        Fal = false(nb,1);
    end
    % Sub一般大于1，所以不会进入第一个 if
    if length(Sub) == 1 && Sub == 0
        % All the cover sets
        ns = nb;
        if nargin == 3
            Sub = true(nb,1);
        else
            Sub = ~Fal;
        end
    % 判断Sub矩阵,是否是布尔类型;;  islogical判断变量是否为逻辑类型
    % 如果Sub不是逻辑类型,则进入if,并将其转换为逻辑类型
    % Sub一般一定是布尔类型，则不进入这个 if
    elseif ~islogical(Sub)
        % Subset of cover sets
        ns = length(Sub);
        if nargin == 3
            sub = false(nb,1);
        else
            sub = Fal;
        end
        sub(Sub) = true;
        Sub = sub;
    else
        % Subset of cover sets
        % 返回Sub中,非零的个数
        ns = nnz(Sub);
    end
    
    % ns是 集合中未被连接的个数，创立一个元胞数组
    Components = cell(ns,1);
    % 根据ns数量，生成一个全0矩阵
    CompSize = zeros(ns,1,'uint32');
    nc = 0;      % number of components found
    m = 1;
    % m 位置的集合 已经被连接的集合
    % 若 ~sub(m)==1, 表示该处集合已连结，若 ~sub(m)==0, 表示该处集合未连结
    while ~Sub(m)
        m = m+1;
    end
    i = 0;
    Comp = zeros(ns,1,'uint32');
    % ns是非零总数
    while i < ns
        % Add 是 m处 的邻居矩阵
        Add = Nei{m};
        % I 是Add在Sub中是什么状态，是否被连结
        I = Sub(Add);
        % Add 保留下 m的邻居矩阵
        Add = Add(I);
        % a 是找到的新邻居个数
        a = length(Add);
        % 在comp第一个位置，保存 m
        Comp(1) = m;
        % 设为0 证明已经被连结
        Sub(m) = false;
        t = 1;
        while a > 0
            % 将Add的内容，添加到m的后面
            Comp(t+1:t+a) = Add;
            % 变为 false 表示已连接
            Sub(Add) = false;
            % 确认当前索引值 comp内的
            t = t+a;
            % 找到Add的所有邻居矩阵（存在重复）
            Add = vertcat(Nei{Add});
            I = Sub(Add);
            Add = Add(I);
            % select the unique elements of Add:
            % n是Add 筛选后的长度
            n = length(Add);
            if n > 2
                % 长度为n全为1的向量
                I = true(n,1);
                for j = 1:n
                    % Fal内的元素理论上均为0，所以第一次一定进入第一个if内
                    % 若Add内有重复的元素，那么I内就会变为false
                    if ~Fal(Add(j))
                        Fal(Add(j)) = true;
                    else
                        I(j) = false;
                    end
                end
                % 将Fal内重置为均为0
                Fal(Add) = false;
                % Add内为 无重复值（所以为什么不使用unique函数）
                Add = Add(I);
            elseif n == 2
                if Add(1) == Add(2)
                    Add = Add(1);
                end
            end
            % a是筛选后Add的长度
            a = length(Add);
        end
        i = i+t;
        % MinSize = 1, 所以一定进入第一个if
        if t >= MinSize
            % nc是新找到的集合数
            nc = nc+1;
            % 将新找到的这些邻居,保存在Components内
            Components{nc} = uint32(Comp(1:t));
            % 保存新邻居的个数
            CompSize(nc) = t;
        end
        if i < ns
            while m <= nb && Sub(m) == false
                m = m+1;
            end
        end
        % 既然每次Comp, Components, ComSize都会重新赋值,那为什么还要初始化设为一个特别大的向量呢
        % 我猜是为了防止,剩下的所有集合,被一次性全部找到
    end
    Components = Components(1:nc);
    CompSize = CompSize(1:nc);
else
    Components = cell(0,1);
    CompSize = 0;
end