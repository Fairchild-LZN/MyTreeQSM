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

function segment = segments(cover,Base,Forb)

% 返回segment结构体
% 其中包括
% segments元组
% ParentSegment向量，保存他们的父节点，每个圆柱体只有一个父亲
% ChildSegment元组，保存他们的孩子节点，每个圆柱体可以有多个孩子

% ---------------------------------------------------------------------
% SEGMENTS.M        Segments the covered point cloud into branches.
%
% Version 2.10
% Latest update     16 Aug 2017
%
% Copyright (C) 2013-2017 Pasi Raumonen
% ---------------------------------------------------------------------

% Segments the tree into branches and records their parent-child-relations. 
% Bifurcations are recognized by studying connectivity of a "study"
% region moving along the tree. In case of multiple connected components 
% in "study", the components are classified as the continuation and branches.
%
% Inputs:
% cover         Cover sets
% Base          Base of the tree
% Forb          Cover sets not part of the tree
%
% Outputs:
% segment       Structure array containing the followin fields:
%   segments          Segments found, (n_seg x 1)-cell, each cell contains a cell array the cover sets
%   ParentSegment     Parent segment of each segment, (n_seg x 1)-vector,
%                       equals to zero if no parent segment
%   ChildSegment      Children segments of each segment, (n_seg x 1)-cell

Nei = cover.neighbor;
% nb保存cover的集合数
nb = size(Nei,1);           % The number of cover sets
% a为点云生成的最多segment的个数
a = max([200000 nb/100]);   % Estimate for maximum number of segments
SBas = cell(a,1);           % The segment bases found
Segs = cell(a,1);           % The segments found
SPar = zeros(a,2,'uint32'); % The parent segment of each segment
SChi = cell(a,1);           % The children segments of each segment

% Initialize SChi
% 初始化Segment的孩子节点5000x1的全0矩阵
SChi{1} = zeros(5000,1,'uint32');
% 其余初始化为200x1的全0矩阵
C = zeros(200,1);
for i = 2:a
    SChi{i} = C;
end
NChi = zeros(a,1);      % Number of child segments found for each segment

Fal = false(nb,1);      % Logical false-vector for cover sets
s = 1;                  % The index of the segment under expansion
% b中保存最后一次开始寻找时的 根的 索引值
% 0329更新，每发现一个新的分支b+1
b = s;                  % The index of the latest found base

SBas{s} = Base;
Seg = cell(1000,1);    % The cover set layers in the current segment
Seg{1} = Base;

ForbAll = Fal;       % The forbidden sets
ForbAll(Forb) = true;
ForbAll(Base) = true;
Forb = ForbAll;      % The forbidden sets for the segment under expansion

Continue = true; % True as long as the component can be segmented further 
NewSeg = true;   % True if the first Cut for the current segment
nl = 1;          % The number of cover set layers currently in the segment

% Segmenting stops when there are no more segments to be found
% 当没有新的seg可以找到时停止
while Continue && (b < nb)
    
    % Update the forbidden sets
    % 将已经拼好的seg剔除
    Forb(Seg{nl}) = true;
    
    % Define the study
    % Nei是邻居集合
    % Seg是
    % Forb是base和不属于树的点集
    % Fal是全0矩阵
    Cut = define_cut(Nei,Seg{nl},Forb,Fal);
    % 邻居集合内的数量
    CutSize = length(Cut);
    
    if NewSeg
        NewSeg = false;
        % 返回CutSize和6内的最小值
        ns = min(CutSize,6);
    end
    
    % Define the components of cut and study regions
    if CutSize > 0
        % Nei是邻居集合
        % Cut新的集合索引
        % CutSize是新集合长度
        % Fal全0
        CutComps = cut_components(Nei,Cut,CutSize,Fal,Fal);
        nc = size(CutComps,1);
        if nc > 1
            % study_components研究哪儿个是主干，哪儿个是分支
            [StudyComps,Bases,CompSize,Cont,BaseSize] = ...
                study_components(Nei,ns,Cut,CutComps,Forb,Fal,Fal);
            nc = length(Cont);
        end
    else
        nc = 0;
    end
    
    % Classify study region components
    if nc == 1
        % One component, continue expansion of the current segment
        % nl代表有多少层，一层一层往上找
        nl = nl+1;
        % 返回Cut的列数
        if size(Cut,2) > 1
            % 如果大于1，则转置
            Seg{nl} = Cut';
        else
            Seg{nl} = Cut;
        end
    elseif nc > 1
        % Classify the components of the Study region
        % 新邻居集合的个数
        Class = component_classification(CompSize,Cont,BaseSize,CutSize);
        
        for i = 1:nc
            % 如果是分支class==1，class==0是主干
            if Class(i) == 1
                % 取出新的分支的Base底部
                Base = Bases{i};
                % ForbAll为true代表，已被记录
                ForbAll(Base) = true;
                Forb(StudyComps{i}) = true;
                % 将cut内属于分支的部分去除掉
                J = Forb(Cut);
                Cut = Cut(~J);
                % 发现新的分支
                % 记录新的base
                b = b+1;
                SBas{b} = Base;
                % 0328 s的意思还不懂
                % 0329 盲猜这个s不会指的是第几个分支？？
                % 0331 s=1代表主干（第一个分支），s=2代表第二个base开始的树枝
                % nl代表它的父节点是在当前分支的第几层
                SPar(b,:) = [s nl];
                % 发现的子分支个数
                NChi(s) = NChi(s)+1;
                % 0331 这步可能是为了后面找下一个base用的
                % s分支（主干也是一种分支），的第Nchi(s)个分支，的b位置的base（上面将base保存在sBas(b)中的）
                SChi{s}(NChi(s)) = b;
            end
        end
        
        % Define the new cut.
        % If the cut is empty, determine the new base
        % 0331 下面是瞎猜的，调试没成功（有问题！）
        % 因为上面的for循环，将cut内属于分支的部分都去除掉
        % 如果（到头了？），这样cut内的所有元素都没了，就有可能进入下面这个if中
        if isempty(Cut)
            Segs{s} = Seg(1:nl);
            S = vertcat(Seg{1:nl});
            ForbAll(S) = true;

            if s < b
                s = s+1;
                Seg{1} = SBas{s};
                Forb = ForbAll;
                NewSeg = true;
                nl = 1;
            else
                Continue = false;
            end
        else
            if size(Cut,2) > 1
                Cut = Cut';
            end
            nl = nl+1;
            Seg{nl} = Cut;
        end
    
    % 当define_cut返回空，即cutsize==0时，会进入这个分支
    else
        % If the study region has zero size, then the current segment is
        % complete and determine the base of the next segment
        % 如果当前的分支没有新的邻居，
        % 那么找到Base开始另一个分支

        % 将这个分支，每一个layer(每层)的集合索引，都保存在segs元胞数组中的一个元素内
        Segs{s} = Seg(1:nl);
        % 取出Segs中的每个元素
        S = vertcat(Seg{1:nl});
        % 标记为true，表示已被标记
        ForbAll(S) = true;
        
        % s是当前正在进行的分支索引，b是现在一共发现了多少个分支个数
        if s < b
            s = s+1;
            % 将新的Base赋值给Seg的底部
            Seg{1} = SBas{s};
            % Forb中为1则表示已被标记
            Forb = ForbAll;
            NewSeg = true;
            nl = 1;
        else
            Continue = false;
        end
    end
end
% 取出每个分支的Segs
Segs = Segs(1:b);
% 取出每个分支的父节点，第一列代表是第几个分支，第二列代表在这个分支的第几层
SPar = SPar(1:b,:);
% schi是元胞数组，里面第n个数组（第n个分支）内保存的是在SPar内的索引
schi = SChi(1:b);

% Define output
SChi = cell(b,1);
for i = 1:b
    if NChi(i) > 0
        % 转换成int类型
        SChi{i} = uint32(schi{i}(1:NChi(i)));
    else
        SChi{i} = zeros(0,1,'uint32');
    end
    S = Segs{i};
    % 将Segs中保存的转换为int类型
    for j = 1:size(S,1)
        S{j} = uint32(S{j});
    end
    Segs{i} = S;
end
clear Segment
segment.segments = Segs;
segment.ParentSegment = SPar;
segment.ChildSegment = SChi;

end % End of the main function


% Define subfunctions

function Cut = define_cut(Nei,CutPre,Forb,Fal)

% Defines the "Cut" region
% Cut保存，当前集合内的全部邻居集合
Cut = vertcat(Nei{CutPre});
Cut = unique_elements(Cut,Fal);
% 删除已经组成seg部分的点
% 每次是从下往上找的，把下面已经分好seg的集合去除掉
I = Forb(Cut);
Cut = Cut(~I);
end % End of function 


function [Components,CompSize] = cut_components(Nei,Cut,CutSize,Fal,False)

% Define the connected components of the Cut
if CutSize == 1
    % Cut is connected and therefore Study is also
    CompSize = 1;
    Components = cell(1,1);
    Components{1} = Cut;
elseif CutSize == 2
    I = Nei{Cut(1)} == Cut(2);
    if any(I)
        Components = cell(1,1);
        Components{1} = Cut;
        CompSize = 1;
    else
        Components = cell(2,1);
        Components{1} = Cut(1);
        Components{2} = Cut(2);
        CompSize = [1 1];
    end
elseif CutSize == 3
    I = Nei{Cut(1)} == Cut(2);
    J = Nei{Cut(1)} == Cut(3);
    K = Nei{Cut(2)} == Cut(3);
    if any(I)+any(J)+any(K) >= 2
        CompSize = 1;
        Components = cell(1,1);
        Components{1} = Cut;
    elseif any(I)
        Components = cell(2,1);
        Components{1} = Cut(1:2);
        Components{2} = Cut(3);
        CompSize = [2 1];
    elseif any(J)
        Components = cell(2,1);
        Components{1} = Cut([1 3]');
        Components{2} = Cut(2);
        CompSize = [2 1];
    elseif any(K)
        Components = cell(2,1);
        Components{1} = Cut(2:3);
        Components{2} = Cut(1);
        CompSize = [2 1];
    else
        CompSize = [1 1 1];
        Components = cell(3,1);
        Components{1} = Cut(1);
        Components{2} = Cut(2);
        Components{3} = Cut(3);
    end
else
    Components = cell(CutSize,1);
    CompSize = zeros(CutSize,1);
    Comp = zeros(CutSize,1);
    Fal(Cut) = true;
    nc = 0;      % number of components found
    m = Cut(1);
    i = 0;
    while i < CutSize
        % 找到当前集合m的所有邻居集合
        % 找到的邻居集合，一定也在Cut内，但会多出已经被seg的集合
        Added = Nei{m};
        % 只保留cut内的集合
        I = Fal(Added);
        Added = Added(I);
        % 计算去除后的集合个数
        a = length(Added);
        Comp(1) = m;
        % 在Fal内标记该集合以被Seg
        Fal(m) = false;
        t = 1;
        % 下面这个while有点类似，沿着圆柱xoy平面横着找邻居集合
        while a > 0
            % 将邻居集合保存在Comp内
            Comp(t+1:t+a) = Added;
            % 在Fal内标记m的邻居集合被seg
            Fal(Added) = false;
            t = t+a;
            Ext = vertcat(Nei{Added});
            Ext = unique_elements(Ext,False);
            I = Fal(Ext);
            Added = Ext(I);
            a = length(Added);
        end
        i = i+t;
        % 如果nc>1则代表碰见分支了
        nc = nc+1;
        Components{nc} = Comp(1:t);
        CompSize(nc) = t;
        if i < CutSize
            J = Fal(Cut);
            m = Cut(J);
            m = m(1);
        end
    end
    Components = Components(1:nc);
    CompSize = CompSize(1:nc);
end

end % End of function


function [Components,Bases,CompSize,Cont,BaseSize] = ...
    study_components(Nei,ns,Cut,CutComps,Forb,Fal,False)
% 输入值
% (ns-1)是算分支数量
% cut是全部邻居集合
% cutcomps是每个分支的邻居集合
% Fal是全0矩阵
% False是全0矩阵

% 输出值
% components 元胞数组类型，保存每个分支，向外延申6圈内的集合索引
% Bases 元胞数组类型，保存每个分支，底部集合索引（Bases所有元素的和 == Cut内元素）
% CompSize int类型，保存每个分支，向外延申6圈内的集合个数
% Cont 逻辑类型，判断分支是否真的为分支，如果是错误分类进去的，则返回false
% BaseSize int类型，保存每个分支，底部集合个数

% Define Study as a cell-array
Study = cell(ns,1);
StudySize = zeros(ns,1);
Study{1} = Cut;
% 总个数
StudySize(1) = length(Cut);
if ns >= 2
    N = Cut;
    i = 1;
    while i < ns
        Forb(N) = true;
        N = vertcat(Nei{N});
        % 去除重复值，去除Fal内为1的值
        N = unique_elements(N,Fal);
        I = Forb(N);
        N = N(~I);
        % 如果N中不为空
        if ~isempty(N)
            i = i+1;
            % 将i的邻居集合保存进来
            Study{i} = N;
            % 保存邻居集合的数量
            StudySize(i) = length(N);
        else
            Study = Study(1:i);
            StudySize = StudySize(1:i);
            i = ns+1;
        end
    end
end

% Define study as a vector
% 统计长度，默认最大为6
ns = length(StudySize);
% 求邻居集合的总数
studysize = sum(StudySize);
% 将所有邻居放在一起
study = vertcat(Study{:});

% Determine the components of study
% 分支+主干个数 = cutcomp
nc = size(CutComps,1);
i = 1; % index of cut component
% j表示有多少个集合被统计记录
j = 0; % number of elements attributed to components
% k表示需要分辨的树干的个数
k = 0; % number of study components

% ***
% 这里就将下面循环限制在study内的集合了
% ***
Fal(study) = true;
Components = cell(nc,1);
CompSize = zeros(nc,1);
Comp = zeros(studysize,1);
while i <= nc
    % 先去除某段“树干”
    C = CutComps{i};
    while j < studysize
        a = length(C);
        % 将C中保存的集合，分别赋值给1-a的位置
        Comp(1:a) = C;
        % false表示已被记录
        Fal(C) = false;
        if a > 1
            % 往外延申一圈
            Add = unique_elements(vertcat(Nei{C}),False);
        else
            Add = Nei{C};
        end
        t = a;
        % 去除已被记录的集合
        I = Fal(Add);
        Add = Add(I);
        % 取新集合长度
        a = length(Add);
        while a > 0
            Comp(t+1:t+a) = Add;
            Fal(Add) = false;
            t = t+a;
            Add = vertcat(Nei{Add});
            Add = unique_elements(Add,False);
            % Fal把这个while的研究集合限制在study内了
            I = Fal(Add);
            Add = Add(I);
            a = length(Add);
        end
        j = j+t;
        k = k+1;
        % 这个分支，包括的集合索引
        Components{k} = Comp(1:t);
        % 这个分支，包括的集合个数
        CompSize(k) = t;
        if j < studysize
            % c被初始化为0矩阵（归零）
            C = zeros(0,1);
            while i < nc && isempty(C)
                i = i+1;
                % 将下一个“树干”赋值给c开启下一段查找
                C = CutComps{i};
                J = Fal(C);
                C = C(J);
            end
            % 如果 i==nc 则所有“树干”查找完毕
            if i == nc && isempty(C)
                j = studysize;
                i = nc+1;
            end
        else
            i = nc+1;
        end
    end
    % 将生成的保存
    % 记录每个“树干”的集合索引
    Components = Components(1:k);
    % 记录每个“树干”内额个数
    CompSize = CompSize(1:k);
end

% Determine BaseSize and Cont
Cont = true(k,1);
BaseSize = zeros(k,1);
Bases = cell(k,1);
if k > 1
    Forb(study) = true;
    Fal(study) = false;
    Fal(Study{1}) = true;
    for i = 1:k
        % Determine the size of the base of the components
        Set = unique_elements([Components{i}; Study{1}],False);
        False(Components{i}) = true;
        % set在False内为1，并且，Set在Fal内为1
        % I是同时属于Components和Study的点
        % Components是分“树干”时保存的集合
        % Study是一开始，进行的六轮外扩的第一轮集合
        % Study是不区分“树干”的往外扩，所以主干和分支都在往外扩
        % 下面这步为的是，去除第一轮往外扩获得的集合中，属于这个分支的集合索引
        I = False(Set)&Fal(Set);
        % 归零
        False(Components{i}) = false;
        Set = Set(I);
        Bases{i} = Set;
        BaseSize(i) = length(Set);
    end
    Fal(Study{1}) = false;
    Fal(Study{ns}) = true;
    Forb(study) = true;
    % 和上面的循环逻辑一样
    % 目的是判断在最外面一圈（第六圈），是否也有相同的集合元素
    for i = 1:k
        % Determine if the component can be extended
        Set = unique_elements([Components{i}; Study{ns}],False);
        False(Components{i}) = true;
        I = False(Set)&Fal(Set);
        False(Components{i}) = false;
        Set = Set(I);
        if ~isempty(Set)
            N = vertcat(Nei{Set});
            N = unique_elements(N,False);
            I = Forb(N);
            N = N(~I);
            if isempty(N)
                Cont(i) = false;
            end
        else
            Cont(i) = false;
        end
    end
end

end % End of function


function Class = component_classification(CompSize,Cont,BaseSize,CutSize)
% 输入量

% Classifies study region components:

% 如果Class == 0 表示是树的主干
% 如果Class == 1 表示是树的分支
% Class(i) == 0 continuation
% Class(i) == 1 branch

% 需要判断的分支个数
nc = size(CompSize,1);
StudySize = sum(CompSize);
% nc x 1的全1矩阵
Class = ones(nc,1);     % true if a component is a branch to be further segmented
ContiComp = 0;
% Simple initial classification
for i = 1:nc
    % 如果基部==总体
    % 自己就是树的主干
    if BaseSize(i) == CompSize(i) && ~Cont(i)
        % component has no expansion, not a branch
        Class(i) = 0;
    % 如果基部有1个集合，总体有2个集合
    % 那么是个很小的外扩，不算一个分支
    elseif BaseSize(i) == 1 && CompSize(i) <= 2 && ~Cont(i)
        % component has very small expansion, not a branch
        Class(i) = 0;
    elseif BaseSize(i)/CutSize < 0.05 && 2*BaseSize(i) >= CompSize(i) && ~Cont(i)
        % component has very small expansion or is very small, not a branch
        Class(i) = 0;
    % 外扩的很少不算一个分支
    elseif CompSize(i) <= 3 && ~Cont(i)
        % very small component, not a branch
        Class(i) = 0;
    % 在Base中占70%，或，在六轮外扩的总数占总体的70%
    elseif BaseSize(i)/CutSize >= 0.7 || CompSize(i) >= 0.7*StudySize
        % continuation of the segment
        Class(i) = 0;
        ContiComp = i;
    % 如果不是以上情况，那么它可能就是一个分支
    else
        % Component is probably a branch
    end
end

% 如果每个分支，经过上述这些分类，都没有被认定为主干，那么就进入下面的if
Branches = Class == 1;
if ContiComp == 0 && any(Branches)
    Ind = (1:1:nc)';
    Branches = Ind(Branches);
    % 选择经过六轮后，包含的集合个数最多的分支，将它定位树的主干
    [~,I] = max(CompSize(Branches));
    Class(Branches(I)) = 0;
end

end % End of function
