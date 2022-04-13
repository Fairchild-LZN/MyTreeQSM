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

function RS = relative_size(P,cover,segment)

% ---------------------------------------------------------------------
% RELATIVE_SIZE.M   Determines relative cover set size for points in new covers
%
% Version 2.00
% Latest update     16 Aug 2017
%
% Copyright (C) 2014-2017 Pasi Raumonen
% ---------------------------------------------------------------------
% 
% Uses existing segmentation and its branching structure to determine
% relative size of the cover sets distributed over new covers. The idea is 
% to decrease the relative size as the branch size decreases. This is 
% realised so that the relative size at the base of a branch is
% proportional to the size of the stem's base, measured as number of
% cover sets in the first few layers. Also when we approach the
% tip of the branch, the relative size decreases to the minimum. 
% Maximum relative size is 256 at the bottom of the
% stem and the minimum is 1 at the tip of every branch.
%
% Output:
% RS    Relative size (1-256), uint8-vector, (n_points x 1)

Bal = cover.ball;
Cen = cover.center;
Nei = cover.neighbor;
Segs = segment.segments;
SChi = segment.ChildSegment;

% 点云xyz的个数
np = size(P,1);     % number of points
% Segs段落的个数
ns = size(Segs,1);  % number of segments

% APRIORI算法???

%% Use branching order and height as apriori info
% Determine the branch orders of the segments
Ord = zeros(ns,1);
C = SChi{1};
order = 0;
while ~isempty(C)
    order = order+order;
    Ord(C) = order;
    C = vertcat(SChi{C});
end
maxO = order+1; % maximum branching order (plus one)

% Determine tree height
% 整棵树最高的点
Top = max(P(Cen,3));
% 整棵树最低的点
Bot = min(P(Cen,3));
% 树的高度
H = Top-Bot;

%% Determine "base size" compared to the stem base
% BaseSize is the relative size of the branch base compared to the stem
% base, measured as number of cover sets in the first layers of the cover
% sets. If it is larger than apriori upper limit based on branching order
% and branch height, then correct to the apriori limit 
BaseSize = zeros(ns,1);
% Determine first the base size at the stem
S = Segs{1};
n = size(S,1);
if n >= 2
    m = min([6 n]);
    BaseSize(1) = mean(cellfun(@length,S(2:m)));
else
    BaseSize(1) = length(S{1});
end
% Then define base size for other segments
for i = 2:ns
    S = Segs{i};
    n = size(S,1);
    if n >= 2
        m = min([6 n]);
        % 根据每个树枝，所对应的每层layer包括的集合个数，去区分分段
        BaseSize(i) = ceil(mean(cellfun(@length,S(2:m)))/BaseSize(1)*256);
    else
        BaseSize(i) = length(S{1})/BaseSize(1)*256;
    end
    % 找到这些集合最小的中心点高度
    bot = min(P(Cen(S{1}),3));
    % 当前的高度 
    h = bot-Bot; % height of the segment's base
    % “反向”与整棵树的比    (H-h)/H
    BS = ceil(256*(maxO-Ord(i))/maxO*(H-h)/H); % maximum apriori base size
    if BaseSize(i) > BS
        BaseSize(i) = BS;
    end
end
BaseSize(1) = 256;

%% Determine relative size for points
TS = 1;
RS = zeros(np,1,'uint8');
for i = 1:ns
    S = Segs{i};
    s = size(S,1);
    for j = 1:s
        Q = S{j};
        RS(vertcat(Bal{Q})) = BaseSize(i)-(BaseSize(i)-TS)*sqrt((j-1)/s);
    end
end

%% Adjust the relative size at the base of child segments
% 备份
RS0 = RS;
for i = 1:ns
    C = SChi{i};
    n = length(C);
    if n > 0
        for j = 1:n
            % 选出当前的子节点分支
            S = Segs{C(j)};
            % 取出第一层Base
            B = S{1};
            % 找邻居集合（上下都有）
            N = vertcat(Nei{B});
            if size(S,1) > 1
                % 返回N中存在，而S{2}中不存在的
                % 差集
                % N之前是，S第一层集合的“邻居”集合
                % 排除掉不属于S第二层的集合
                % 剩下的可能是属于父节点的吧？0408
                N = setdiff(N,S{2});
            end
            % 返回N和B的并集
            % S{1}和属于父节点的集合
            N = union(N,B);
            N = vertcat(Bal{N});
            % 0408 不理解啊，为啥base要除以2
            % 0412
            % 前面计算时，是从2：m的，就代表没有用第一层
            % 再结合correct_segment内的最后一步
            % 将每个分支都“往后延”了一截（这步仍然存疑）
            % 所以这里骚操作了一波（但为什么除2仍不理解）
            % 0412再更新
            % 上面好像又不对
            % 主要是setdiff和union这两步到底去除添加了哪儿些集合？
            RS(N) = RS0(N)/2;
        end
    end
end
