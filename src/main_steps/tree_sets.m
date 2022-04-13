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

function [cover,Base,Forb] = tree_sets(P,cover,inputs,segment)

% ---------------------------------------------------------------------
% TREE_SETS.M       Determines the base of the trunk and the cover sets 
%                   belonging to the tree, updates the neighbor-relation
%
% Version 2.2.0
% Latest update     13 Aug 2020
%
% Copyright (C) 2013-2020 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Determines the cover sets that belong to the tree. Determines also the
% base of the tree and updates the neighbor-relation such that all of the
% tree is connected, i.e., the cover sets belonging to the tree form a
% single connected component. Optionally uses information from existing
% segmentation to make sure that stem and 1st-, 2nd-, 3rd-order branches 
% are properly connnected.
% ---------------------------------------------------------------------
% Inputs:
% P             Point cloud
% cover         Cover sets, their centers and neighbors
% PatchDiam     Minimum diameter of the cover sets
% OnlyTree      Logical value indicating if the point cloud contains only
%                   points from the tree to be modelled
% segment       Previous segments
%
% Outputs:
% cover     Cover sets with updated neigbors
% Base      Base of the trunk (the cover sets forming the base)
% Forb      Cover sets not part of the tree
% ---------------------------------------------------------------------

% Changes from version 2.1.0 to 2.2.0, 13 Aug 2020:  
% 1) "define_base_forb": Changed the base height specification from
%     0.1*aux.Height to 0.02*aux.Height
% 2) "define_base_forb": changed the cylinder fitting syntax corresponding
%     to the new input and outputs of "least_squares_cylinder" 
% 3) "make_tree_connected”: Removed "Trunk(Base) = false;" at the beginning 
%     of the function as unnecessary and to prevent errors in a special case 
%     where the Trunk is equal to Base.
%	4) "make_tree_connected”: Removed from the end the generation of "Trunk" 
%     again and the new call for the function 
%	5) "make_tree_connected”: Increased the minimum distance of a component 
%     to be removed from 8m to 12m.

% Changes from version 2.0.0 to 2.1.0, 11 Oct 2019:  
% 1) "define_main_branches": modified the size of neighborhood "balls0", 
%    added seven lines of code, prevents possible error of too low or big 
%    indexes on "Par"
% 2) Increased the maximum base height from 0.5m to 1.5m
% 3) "make_tree_connected": added at the end a call for the function itself,
%    if the tree is not yet connected, thus running the function again if 
%    necessary 

%% Define auxiliar object
clear aux
% cover_sets的个数;; max返回最大值,size返回矩阵是x*x维的
aux.nb = max(size(cover.center));   % number of cover sets
% aux.nb x 1 全为0的矩阵
aux.Fal = false(aux.nb,1);
% Index, 1到nb,间隔为1的序号
aux.Ind = (1:1:aux.nb)';
% 取出cover.center的中心点坐标
aux.Ce = P(cover.center,1:3); % Coordinates of the center points
% 高度最小值
aux.Hmin = min(aux.Ce(:,3));
% 高度最大值-最小值
aux.Height = max(aux.Ce(:,3))-aux.Hmin;

%% Define the base of the trunk and the forbidden sets
if nargin == 3
    [Base,Forb,cover] = define_base_forb(P,cover,aux,inputs);
else
    inputs.OnlyTree = true;
    [Base,Forb,cover] = define_base_forb(P,cover,aux,inputs,segment);
end

%% Define the trunk (and the main branches)
if nargin == 3
    [Trunk,cover] = define_trunk(cover,aux,Base,Forb,inputs);
else
    [Trunk,cover] = define_main_branches(cover,segment,aux,inputs);
end

%% Update neighbor-relation to make the whole tree connected
[cover,Forb] = make_tree_connected(cover,aux,Forb,Base,Trunk,inputs);

end % End of the main function


function [Base,Forb,cover] = define_base_forb(P,cover,aux,inputs,segment)

% 定义树的底部和不属于树的集合(集合内含有不包括树的点)

% Defines the base of the stem and the forbidden sets (the sets containing
% points not from the tree, i.e, ground, understory, etc.)
Ce = aux.Ce;

% inputs.onlytree(传入的点云数据只有树的点云,没有其他不属于树的杂质)
if inputs.OnlyTree && nargin == 4
    % No ground in the point cloud, the base is the lowest part
    % 定义树的基部,(1.5 or 0.02倍树高的中心点)
    BaseHeight = min(1.5,0.02*aux.Height);

    % 0323
    % **下面这行,是运行前瞎猜的,属实小丑
    % 理论上 I 内均为0; 因为所有中心点肯定大于,最小树高+Base底部
    % **

    % 我不理解为什么会存在中心点,比Hmin还低

    % 0327更新，我理解了，我之前不理解属实犯蠢

    % I中是中心点高度坐标，小于最小值+Base高度
    I = Ce(:,3) < aux.Hmin+BaseHeight;
    % 取出比小于基部所对应的点
    Base = aux.Ind(I);
    Forb = aux.Fal;
    % Make sure the base, as the bottom of point cloud, is not in multiple parts
    % Base即为底部的点坐标

    % ******
    % base中分别取xy轴的最大最小值,并相减,然后取出较大值
    Wb = max(max(Ce(Base,1:2))-min(Ce(Base,1:2)));
    % 中心点中分别取xy轴的最大最小值,并相减,然后取出较大值
    Wt = max(max(Ce(:,1:2))-min(Ce(:,1:2)));
    % ******

    k = 1;
    while k <= 5 && Wb > 0.3*Wt
        BaseHeight = BaseHeight-0.05;
        BaseHeight = max(BaseHeight,0.05);
        if BaseHeight > 0
            I = Ce(:,3) < aux.Hmin+BaseHeight;
        else
            [~,I] = min(Ce(:,3));
        end
        Base = aux.Ind(I);
        Wb = max(max(Ce(Base,1:2))-min(Ce(Base,1:2)));
        k = k+1;
    end
elseif inputs.OnlyTree
    % Select the stem sets from the previous segmentation and define the
    % base
    BaseHeight = min(1.5,0.02*aux.Height);

    % SegmentOfPoint对应的是P中每个点所对应的Seg索引值
    SoP = segment.SegmentOfPoint(cover.center);
    % 取出属于主干的点
    stem = aux.Ind(SoP == 1);
    % 取出小于底部base高度的值
    I = Ce(stem,3) < aux.Hmin+BaseHeight;
    Base = stem(I);
    % 全0
    Forb = aux.Fal;
else
    % Point cloud contains non-tree points.
    % Determine the base from the "height" and "density" of cover sets
    % by projecting the sets to the xy-plane
    Bal = cover.ball;
    Nei = cover.neighbor;
    
    % The vertices of the rectangle containing C
    Min = double(min(Ce));
    Max = double(max(Ce(:,1:2)));
    
    % Number of rectangles with edge length "E" in the plane
    E = min(0.1,0.015*aux.Height);
    n = double(ceil((Max(1:2)-Min(1:2))/E)+1);
    
    % Calculates the rectangular-coordinates of the points
    px = floor((Ce(:,1)-Min(1))/E)+1;
    py = floor((Ce(:,2)-Min(2))/E)+1;
    
    % Sorts the points according a lexicographical order
    LexOrd = [px py-1]*[1 n(1)]';
    [LexOrd,SortOrd] = sort(LexOrd);
    
    Partition = cell(n(1),n(2));
    hei = zeros(n(1),n(2)); % "height" of the cover sets in the squares
    den = hei;  % density of the cover sets in the squares
    baseden = hei;
    p = 1; % The index of the point under comparison
    while p <= aux.nb
        t = 1;
        while (p+t <= aux.nb) && (LexOrd(p) == LexOrd(p+t))
            t = t+1;
        end
        q = SortOrd(p);
        J = SortOrd(p:p+t-1);
        Partition{px(q),py(q)} = J;
        p = p+t;
        K = ceil(10*(Ce(J,3)-Min(3)+0.01)/(aux.Height-0.01));
        B = K <= 2;
        K = unique(K);
        hei(px(q),py(q)) = length(K)/10;
        den(px(q),py(q)) = t;
        baseden(px(q),py(q)) = nnz(B);
    end
    den = den/max(max(den));  % normalize
    baseden = baseden/max(max(baseden));
    
    % function whose maximum determines location of the trunk
    f = den.*hei.*baseden;  
    % smooth the function by averaging over 8-neighbors
    x = zeros(n(1),n(2));
    y = zeros(n(1),n(2));
    for i = 2:n(1)-1
        for j = 2:n(2)-1
            f(i,j) = mean(mean(f(i-1:i+1,j-1:j+1)));
            x(i,j) = Min(1)+i*E;
            y(i,j) = Min(2)+j*E;
        end
    end
    f = f/max(max(f));
    
    % Trunk location is around the maximum f-value
    I = f > 0.5;
    Trunk0 = Partition(I); % squares that contain the trunk
    Trunk0 = vertcat(Trunk0{:});
    HBottom = min(Ce(Trunk0,3));
    I = Ce(Trunk0,3) > HBottom+min(0.02*aux.Height,0.3);
    J = Ce(Trunk0,3) < HBottom+min(0.08*aux.Height,1.5);
    I = I&J; % slice close to bottom should contain the trunk
    Trunk = Trunk0(I);
    Trunk = union(Trunk,vertcat(Nei{Trunk})); % Expand with neighbors
    Trunk = union(Trunk,vertcat(Nei{Trunk})); % Expand with neighbors
    Trunk = union(Trunk,vertcat(Nei{Trunk})); % Expand with neighbors
    
    % Define connected components of Trunk and select the largest component
    [Comp,CS] = connected_components(Nei,Trunk,0,aux.Fal);
    [~,I] = max(CS);
    Trunk = Comp{I};
    
    % Fit cylinder to Trunk
    I = Ce(Trunk,3) < HBottom+min(0.1*aux.Height,2); % Select the bottom part
    Trunk = Trunk(I);
    Trunk = union(Trunk,vertcat(Nei{Trunk}));
    Points = Ce(Trunk,:);
    c.start = mean(Points);
    c.axis = [0 0 1];
    c.radius = mean(distances_to_line(Points,c.axis,c.start));
    c = least_squares_cylinder(Points,c);
    
    % Remove far away points and fit new cylinder
    dis = distances_to_line(Points,c.axis,c.start);
    [~,I] = sort(abs(dis));
    I = I(1:ceil(0.9*length(I)));
    Points = Points(I,:);
    Trunk = Trunk(I);
    c = least_squares_cylinder(Points,c);
    
    % Select the sets in the bottom part of the trunk and remove sets too
    % far away form the cylinder axis (also remove far away points from sets)
    I = Ce(Trunk0,3) < HBottom+min(0.04*aux.Height,0.6);
    TrunkBot = Trunk0(I);
    TrunkBot = union(TrunkBot,vertcat(Nei{TrunkBot}));
    TrunkBot = union(TrunkBot,vertcat(Nei{TrunkBot}));
    n = length(TrunkBot);
    Keep = true(n,1); % Keep sets that are close enough the axis
    a = max(0.06,0.2*c.radius);
    b = max(0.04,0.15*c.radius);
    for i = 1:n
        d = distances_to_line(Ce(TrunkBot(i),:),c.axis,c.start);
        if d < c.radius+a
            B = Bal{Trunk(i)};
            d = distances_to_line(P(B,:),c.axis,c.start);
            I = d < c.radius+b;
            Bal{Trunk(i)} = B(I);
        else
            Keep(i) = false;
        end
    end
    TrunkBot = TrunkBot(Keep);
    
    % Select the above part of the trunk and combine with the bottom
    I = Ce(Trunk0,3) > HBottom+min(0.03*aux.Height,0.45);
    Trunk = Trunk0(I);
    Trunk = union(Trunk,vertcat(Nei{Trunk}));
    Trunk = union(Trunk,TrunkBot);
    
    BaseHeight = min(1.5,0.02*aux.Height);
    % Determine the base
    Bot = min(Ce(Trunk,3));
    J = Ce(Trunk,3) < Bot+BaseHeight;
    Base = Trunk(J);
    
    % Determine "Forb", i.e, ground and non-tree sets by expanding Trunk
    % as much as possible
    Trunk = union(Trunk,vertcat(Nei{Trunk}));
    Forb = aux.Fal;
    Ground = setdiff(vertcat(Nei{Base}),Trunk);
    Ground = setdiff(union(Ground,vertcat(Nei{Ground})),Trunk);
    Forb(Ground) = true;
    Forb(Base) = false;
    Add = Forb;
    while any(Add)
        Add(vertcat(Nei{Add})) = true;
        Add(Forb) = false;
        Add(Trunk) = false;
        Forb(Add) = true;
    end
    
    % Try to expand the "Forb" more by adding all the bottom sets
    Bot = min(Ce(Trunk,3));
    Ground = Ce(:,3) < Bot+0.03*aux.Height;
    Forb(Ground) = true;
    Forb(Trunk) = false;
    cover.ball = Bal;
end

end % End of function


function [Trunk,cover] = define_trunk(cover,aux,Base,Forb,inputs)
% 0323
% 从Base开始，往上找树干，定义树干为base往上0.25m内的集合均算

% 0327更新
% 从Base开始，往上找树干，每次寻找的高度为0.25m
% 找主干，不在意分支？ 因为有一步要求找到的集合要大于 H 值

% cover 聚类集合
% aux 参数集
% Base 底部集合
% Forb 全0矩阵
% inputs 输入参数

% This function tries to make sure that likely "route" of the trunk from
% the bottom to the top is connected. However, this does not mean that the
% final trunk follows this "route".

Nei = cover.neighbor;
Ce = aux.Ce;
% Determine the output "Trunk" which indicates which sets are part of
% likely trunk

% Trunk内初始化全0
Trunk = aux.Fal;
% 位于Base的集合定义为1
Trunk(Base) = true;
% Expand Trunk from the base above with neighbors as long as possible
Exp = Base; % the current "top" of Trunk
% select the unique neighbors of Exp

% 删掉Nei(Exp)中重复的集合序号
Exp = unique_elements([Exp; vertcat(Nei{Exp})],aux.Fal);
% I 中 取出Trunk中exp的集合
I = Trunk(Exp);
J = Forb(Exp);
Exp = Exp(~I|~J); % Only non forbidden sets that are not already in Trunk
Trunk(Exp) = true; % Add the expansion Exp to Trunk


% 最大高度差 0.25m
L = 0.25; % maximum height difference in Exp from its top to bottom


% H 是底部这些集合中,最高的高度 - L
% 所有属于base底部的最低的高度值;;;(比这个高度还低的可能就是不属于树的点？)
H = max(Ce(Trunk,3))-L; % the minimum bottom heigth for the current Exp
% true as long as the expansion is possible with original neighbors:
FirstMod = true;
% 只要Exp中不为空
while ~isempty(Exp)
    % Expand Trunk similarly as above as long as possible
    H0 = H;
    Exp0 = Exp;
    % 将Exp与Exp的邻居找到,并放到同一个矩阵内,删掉重复值
    % 我不理解为什么Exp定义这么多次, 差距在那儿里呢
    Exp = union(Exp,vertcat(Nei{Exp}));
    % I 中等于0的就是新找到的集合
    I = Trunk(Exp);
    % 找到Exp中,之前没有选入的,即 新expansion出来的
    Exp = Exp(~I);
    % 选出大于等于 H 的值
    % 如果比 H 小，那么就不属于base内--0327存疑
    I = Ce(Exp,3) >= H;
    % 将属于Base的集合，放在Exp内
    Exp = Exp(I);
    Trunk(Exp) = true;
    if ~isempty(Exp)
        H = max(Ce(Exp,3))-L;
    end
    
    % If the expansion Exp is empty and the top of the tree is still over 5
    % meters higher, then search new neighbors from above

    % "||"和"&&" 是短运算符,
    % 若(isempty(Exp)为1,则不计算后续, (A中没有元素，结果为1，A中有元素，结果为0)
    % 若H < H0+inputs.PatchDiam1/2)为false则不计算后续

    % 0327貌似很难进入下面的if语句--存疑
    if (isempty(Exp) || H < H0+inputs.PatchDiam1/2) && H < aux.Height-5
        
        % Generate rectangular partition of the sets
        if FirstMod
            FirstMod = false;
            % The vertices of the rectangle containing C
            Min = double(min(Ce(:,1:2)));
            Max = double(max(Ce(:,1:2)));
            nb = size(Ce,1);
            
            % Number of rectangles with edge length "E" in the plane
            EdgeLenth = 0.2;
            NRect = double(ceil((Max-Min)/EdgeLenth)+1);
            
            % Calculates the rectangular-coordinates of the points
            px = floor((Ce(:,1)-Min(1))/EdgeLenth)+1;
            py = floor((Ce(:,2)-Min(2))/EdgeLenth)+1;
            
            % Sorts the points according a lexicographical order
            LexOrd = [px py-1]*[1 NRect(1)]';
            [LexOrd,SortOrd] = sort(LexOrd);
            
            Partition = cell(NRect(1),NRect(2));
            p = 1; % The index of the point under comparison
            while p <= nb
                t = 1;
                while (p+t <= nb) && (LexOrd(p) == LexOrd(p+t))
                    t = t+1;
                end
                q = SortOrd(p);
                J = SortOrd(p:p+t-1);
                Partition{px(q),py(q)} = J;
                p = p+t;
            end
        end
        
        % Select the region that is connected to a set above it
        if ~isempty(Exp)
            Region = Exp;
        else
            Region = Exp0;
        end
        
        % Select the minimum and maximum rectangular coordinate of the
        % region
        X1 = min(px(Region));
        if X1 <= 2
            X1 = 3;
        end
        X2 = max(px(Region));
        if X2 >= NRect(1)-1
            X2 = NRect(1)-2;
        end
        Y1 = min(py(Region));
        if Y1 <= 2
            Y1 = 3;
        end
        Y2 = max(py(Region));
        if Y2 >= NRect(2)-1
            Y2 = NRect(2)-2;
        end
        
        % Select the sets in the 2 meter layer above the region
        sets = Partition(X1-2:X2+2,Y1-2:Y2+2);
        sets = vertcat(sets{:});
        K = aux.Fal;
        K(sets) = true; % the potential sets
        I = Ce(:,3) > H;
        J = Ce(:,3) < H+2;
        I = I&J&K;
        I(Trunk) = false; % Must be non-Trunk sets
        SetsAbove = aux.Ind(I);
        
        % Search the closest connection between Region and SetsAbove that
        % is enough upward sloping (angle to the vertical has cosine larger 
        % than 0.7)
        if ~isempty(SetsAbove)
            % Compute the distances and cosines of the connections
            n = length(Region);
            m = length(SetsAbove);
            Dist = zeros(n,m);
            Cos = zeros(n,m);
            for i = 1:n
                V = mat_vec_subtraction(Ce(SetsAbove,:),Ce(Region(i),:));
                Len = sum(V.*V,2);
                v = normalize(V);
                Dist(i,:) = Len';
                Cos(i,:) = v(:,3)';
            end
            I = Cos > 0.7; % select those connection with large enough cosines
            % if not any, search with smaller cosines
            t = 0;
            while ~any(I)
                t = t+1;
                I = Cos > 0.7-t*0.05;
            end
            % Search the minimum distance
            Dist(~I) = 3;
            if n > 1 && m > 1
                [d,I] = min(Dist);
                [~,J] = min(d);
                I = I(J);
            elseif n == 1 && m > 1
                [~,J] = min(Dist);
                I = 1;
            elseif m == 1 && n < 1
                [~,I] = min(Dist);
                J = 1;
            else
                I = 1; % the set in component to be connected
                J = 1; % the set in "trunk" to be connected
            end
            
            % Join to "SetsAbove"
            I = Region(I);
            J = SetsAbove(J);
            % make the connection
            Nei{I} = [Nei{I}; J];
            Nei{J} = [Nei{J}; I];
            
            % Expand "Trunk" again
            Exp = union(Region,vertcat(Nei{Region}));
            I = Trunk(Exp);
            Exp = Exp(~I);
            I = Ce(Exp,3) >= H;
            Exp = Exp(I);
            Trunk(Exp) = true;
            H = max(Ce(Exp,3))-L;
        end
    end
end
% 下面这步，貌似是针对进入上面的小if语句使用的
cover.neighbor = Nei;

end % End of function


function [Trunk,cover] = define_main_branches(cover,segment,aux,inputs)

% If previous segmentation exists, then use it to make the sets in its main
% branches (stem and first (second or even up to third) order branches) 
% connected. This ensures that similar branching structure as in the 
% existing segmentation is possible.

Bal = cover.ball;
Nei = cover.neighbor;
Ce = aux.Ce;
% Determine sets in the main branches of previous segmentation
nb = size(Bal,1);
MainBranches = zeros(nb,1);
SegmentOfPoint = segment.SegmentOfPoint;
% Determine which branch indexes define the main branches
% Seg分支总数
MainBranchIndexes = false(max(SegmentOfPoint),1);
% 分支123和主干,表示已分类
MainBranchIndexes(1) = true;
MainBranchIndexes(segment.branch1indexes) = true;
MainBranchIndexes(segment.branch2indexes) = true;
MainBranchIndexes(segment.branch3indexes) = true;
for i = 1:nb
    BranchInd = nonzeros(SegmentOfPoint(Bal{i}));
    if ~isempty(BranchInd)
        ind = min(BranchInd);
        % 如果属于123分支，那么属于MainBranch
        if MainBranchIndexes(ind)
            MainBranches(i) = min(BranchInd);
        end
    end
end

% Define the trunk sets
Trunk = aux.Fal;
% 标记已经被标记的分支
Trunk(MainBranches > 0) = true;

% Update the neighbors to make the main branches connected
% 根据之前定义的集合的中心点，正方体的边长x3
[Par,CC] = cubical_partition(Ce,3*inputs.PatchDiam2Max,10);
Sets = zeros(aux.nb,1,'uint32');
BI = max(MainBranches);
N = size(Par);
for i = 1:BI
    if MainBranchIndexes(i)
        Branch = MainBranches == i; % The sets forming branch "i"
        % the connected components of "Branch":
        Comps = connected_components(Nei,Branch,1,aux.Fal); 
        n = size(Comps,1);
        % Connect the components to each other as long as there are more than
        % one component
        while n > 1
            for j = 1:n
                comp = Comps{j};
                NC = length(comp);
                
                % Determine branch sets closest to the component
                c = unique(CC(comp,:),'rows');
                m = size(c,1);
                t = 0;
                NearSets = zeros(0,1);
                while isempty(NearSets)
                    NearSets = aux.Fal;
                    t = t+1;
                    for k = 1:m
                        x1 = max(1,c(k,1)-t);
                        x2 = min(c(k,1)+t,N(1));
                        y1 = max(1,c(k,2)-t);
                        y2 = min(c(k,2)+t,N(2));
                        z1 = max(1,c(k,3)-t);
                        z2 = min(c(k,3)+t,N(3));
                        balls0 = Par(x1:x2,y1:y2,z1:z2);
                        if t == 1
                            balls = vertcat(balls0{:});
                        else
                            S = cellfun('length',balls0);
                            I = S > 0;
                            S = S(I);
                            balls0 = balls0(I);
                            stop = cumsum(S);
                            start = [0; stop]+1;
                            for l = 1:length(stop)
                                Sets(start(l):stop(l)) = balls0{l};
                            end
                            balls = Sets(1:stop(l));
                        end
                        I = Branch(balls);
                        balls = balls(I);
                        NearSets(balls) = true;
                    end
                    NearSets(comp) = false; % Only the non-component cover sets
                    NearSets = aux.Ind(NearSets);
                end
                
                % Determine the closest sets for "comp"
                if ~isempty(NearSets)
                    d = pdist2(Ce(comp,:),Ce(NearSets,:));
                    if NC == 1 && length(NearSets) == 1
                        IU = 1; % the set in component to be connected
                        JU = 1; % the set in "trunk" to be connected
                    elseif NC == 1
                        [du,JU] = min(d);
                        IU = 1;
                    elseif length(NearSets) == 1
                        [du,IU] = min(d);
                        JU = 1;
                    else
                        [d,IU] = min(d);
                        [du,JU] = min(d);
                        IU = IU(JU);
                    end
                    
                    % Join to the closest component
                    I = comp(IU);
                    J = NearSets(JU);
                    % make the connection
                    Nei{I} = [Nei{I}; J];
                    Nei{J} = [Nei{J}; I];
                end
            end
            
            Comps = connected_components(Nei,Branch,1,aux.Fal);
            n = size(Comps,1);
        end
    end
end

% Update the neigbors to connect 1st-order branches to the stem
Stem = MainBranches == 1;
Stem = aux.Ind(Stem);
MainBranchIndexes = false(max(SegmentOfPoint),1);
MainBranchIndexes(segment.branch1indexes) = true;
BI = max(segment.branch1indexes);
if isempty(BI)
    BI = 0;
end
for i = 2:BI
    if MainBranchIndexes(i)
        Branch = MainBranches == i;
        Branch = aux.Ind(Branch);
        if ~isempty(Branch)
            Neigbors = MainBranches(vertcat(Nei{Branch})) == 1;
            if ~any(Neigbors)
                d = pdist2(Ce(Branch,:),Ce(Stem,:));
                if length(Branch) > 1 && length(Stem) > 1
                    [d,I] = min(d);
                    [d,J] = min(d);
                    I = I(J);
                elseif length(Branch) == 1 && length(Stem) > 1
                    [d,J] = min(d);
                    I = 1;
                elseif length(Stem) == 1 && length(Branch) > 1
                    [d,I] = min(d);
                    J = 1;
                elseif length(Branch) == 1 && length(Stem) == 1
                    I = 1; % the set in component to be connected
                    J = 1; % the set in "trunk" to be connected
                end
                
                % Join the Branch to Stem
                I = Branch(I);
                J = Stem(J);
                Nei{I} = [Nei{I}; J];
                Nei{J} = [Nei{J}; I];
            end
        end
    end
end
cover.neighbor = Nei;
    
end % End of function


function [cover,Forb] = make_tree_connected(cover,aux,Forb,Base,Trunk,inputs)

% cover是cover_sets内生成的集合信息
% aux是tree_sets里的初始化的一些变量集合
% Forb是全为0的矩阵
% Bass是树的基部集合
% Trunk是树的主干，其中为0的代表该集合是可以和base连通的，为1的代表还没和base连通
% inputs是主程序的输入变量集合

% 让整棵树连通起来

% Update neighbor-relation for whole tree such that the whole tree is one
% connected component

% Nei是每个集合的邻居集合
Nei = cover.neighbor;
% Ce是集合的中心点
Ce = aux.Ce;
% Expand trunk as much as possible
Trunk(Forb) = false;
% 将Trunk赋值给Exp，布尔类型
Exp = Trunk;
% any函数, Exp全为0时返回0
% 根据之前的Trunk内集合，往外延搜索
% 经过while后，仍然有很多集合没有连接上，因为点云的断层/遮掩等原因

% 0327下面while，我觉得目的是延申branch
while any(Exp)
    % 寻找当前Exp的邻居集合，并重新赋值给Exp，赋值为true
    Exp(vertcat(Nei{Exp})) = true;
    % 排除Exp内属于Trunk的集合
    Exp(Trunk) = false;
    % 排除Exp内属于Forb的集合（Forb内保存的内容存疑）
    Exp(Forb) = false;
    % 排除Exp内属于base的集合
    Exp(Base) = false;
    % 将剩余的集合序号，在Trunk内赋值，属于Trunk
    Trunk(Exp) = true;
end

% Define "Other", sets not yet connected to trunk or Forb
% 给other赋值全为1的矩阵，aux.Fal全为0
Other = ~aux.Fal;
Other(Forb) = false;
Other(Trunk) = false;
Other(Base) = false;

% Determine parameters on the extent of the "Nearby Space" and acceptable
% component size
% cell size for "Nearby Space" = k0 times PatchDiam:

% ceil函数,向上取整
% ceil内为3（单位不清楚 m ??）,  PatchDiam1 = 0.08
k0 = min(10,ceil(0.2/inputs.PatchDiam1)); 
% 每当当前集合没有新的邻居时，就扩大搜索范围，每次向外扩半径k0

% current cell size, increases by k0 every time when new connections cannot
% be made:
k = k0; 
if inputs.OnlyTree
    Cmin = 0;
else
    % Cmin是向外扩的最小接收半径（若小于Cmin那么就认为这部分是鬼点）
    % 如果程序属于的点云，含有鬼点，则选择 Cmin = 2（单位不清楚）
    Cmin = ceil(0.1/inputs.PatchDiam1);  % minimum accepted component size, 
    % smaller ones are added to Forb, the size triples every round
end

% Determine the components of "Other"
if any(Other)
    Comps = connected_components(Nei,Other,1,aux.Fal);
    % nc是未连结的集合数
    nc = size(Comps,1);
    NonClassified = true(nc,1);
    %plot_segs(P,Comps,6,1,cover.ball)
    %pause
else
    NonClassified = false;
end

% bottom是 base的最低点
bottom = min(Ce(Base,3));
% repeat search and connecting as long as "Other" sets exists
% 若 NonClassified 中全为0，则返回0
while any(NonClassified) 
  npre = nnz(NonClassified); % number of "Other" sets before new connections
  again = true; % check connections again with same "distance" if true
  
  % Partition the centers of the cover sets into cubes with size k*dmin
  [Par,CC] = cubical_partition(Ce,k*inputs.PatchDiam1);
  % 新的元胞数组
  Neighbors = cell(nc,1);
  % 生成 nc x 2 的全0矩阵  
  Sizes = zeros(nc,2);
  % 生成 nc x 1 的全1矩阵   
  Pass = true(nc,1);
  first_round = true;
  while again
    % Check each component: part of "Tree" or "Forb"
    for i = 1:nc
      % 如果NonClassified 和 Pass 中 同时为1
      if NonClassified(i) && Pass(i)
        comp = Comps{i}; % candidate component for joining to the tree
        
        % If the component is neighbor of forbidden sets, remove it
        % 如果J属于Forb，Forb是不属于树的集合
        J = Forb(vertcat(Nei{comp}));
        % 一般不会进这里
        if any(J)
          NonClassified(i) = false;
          Forb(comp) = true;
          Other(comp) = false;
        else
          % Other wise check nearest sets for a connection
          NC = length(comp);
          if first_round
              
            % Select the cover sets the nearest to the component
            % 找到comp内集合 对应的新的分类集合内的位置，“row”的意思不懂
            c = unique(CC(comp,:),'rows');
            m = size(c,1);
            B = cell(m,1);
            % 针对新的范围找新邻居集合
            for j = 1:m
                balls = Par(c(j,1)-1:c(j,1)+1,...
                    c(j,2)-1:c(j,2)+1,c(j,3)-1:c(j,3)+1);
                B{j} = vertcat(balls{:});
            end
            % 新的邻居集合放在一起
            NearSets = vertcat(B{:});
            % Only the non-component cover sets
            % 本身所属的小集合为true
            aux.Fal(comp) = true;
            % I 中为1的是原comp内小集合，其余新的邻居集合为0
            I = aux.Fal(NearSets);
            % 只保留新的邻居集合
            NearSets = NearSets(~I);
            % Fal复原
            aux.Fal(comp) = false;
            % 去掉重复点
            NearSets = unique(NearSets);
            % 保存
            Neighbors{i} = NearSets;
            if isempty(NearSets)
                Pass(i) = false;
            end
            % No "Other" sets
            % 观察新的邻居集合，是否已经确认连结树干
            I = Other(NearSets);
            % 保留，已经连结的集合
            NearSets = NearSets(~I);
          else
            NearSets = Neighbors{i};
            % No "Other" sets
            I = Other(NearSets);
            NearSets = NearSets(~I);
          end
          
          % Select different class from NearSets
          % I 中表示 是否属于Trunk  
          I = Trunk(NearSets);
          % J 中表示 是否不属于树
          J = Forb(NearSets);
          % 将属于树干的部分放入trunk中
          trunk = NearSets(I); % "Trunk" sets
          forb = NearSets(J); % "Forb" sets
          % 这里的sizes是，未连接的集合数
          if length(trunk) ~= Sizes(i,1) || length(forb) ~= Sizes(i,2)
            % 将当前trunk和forb的个数保存进去
            Sizes(i,:) = [length(trunk) length(forb)];
            
            % If large component is tall and close to ground, then
            % search the connection near the component's bottom
            % NC是之前未连结的集合中，延申找出来的集合合集个数
            if NC > 100
              hmin = min(Ce(comp,3));
              H = max(Ce(comp,3))-hmin;
              if H > 5 && hmin < bottom+5
                I = Ce(NearSets,3) < hmin+0.5;
                NearSets = NearSets(I);
                I = Trunk(NearSets);
                J = Forb(NearSets);
                trunk = NearSets(I); % "Trunk" sets
                forb = NearSets(J); % "Forb" sets
              end
            end
            
            % Determine the closest sets for "trunk"
            if ~isempty(trunk)
                % 计算comp和trunk 在ce中所指代的坐标的两点之间距离
                % d 是 comp x trunk 的二维矩阵，保存的是欧式距离
                d = pdist2(Ce(comp,:),Ce(trunk,:));
                if NC == 1 && length(trunk) == 1
                  dt = d; % the minimum distance
                  IC = 1; % the set in component to be connected
                  IT = 1; % the set in "trunk" to be connected
                elseif NC == 1
                  [dt,IT] = min(d);
                  IC = 1;
                elseif length(trunk) == 1
                  [dt,IC] = min(d);
                  IT = 1;
                else
                  % d 保存之前d中每列的最小值， IC 保存之前d中每列最小值的行号
                  % 找到最小值的行列位置
                  [d,IC] = min(d);
                  [dt,IT] = min(d);
                  IC = IC(IT);
                end
            else
              dt = 700;
            end
            
            % Determine the closest sets for "forb"
            if ~isempty(forb)
              d = pdist2(Ce(comp,:),Ce(forb,:));
              df = min(d);
              if length(df) > 1
                  df = min(df);
              end
            else
              df = 1000;
            end
            
            % Determine what to do with the component
            if (dt > 12 && dt < 100) || (NC < Cmin && dt > 0.5 && dt < 10)
              % 移走很小的component
              % Remove small isolated component
              Forb(comp) = true;
              Other(comp) = false;
              NonClassified(i) = false;
            elseif 3*df < dt || (df < dt && df > 0.25)
              % 将这部分加入到 不属于树的集合内
              % Join the component to "Forb"
              Forb(comp) = true;
              Other(comp) = false;
                NonClassified(i) = false;
            % 如果trunk和forb内均为空
            elseif (df == 1000 && dt == 700) || dt > k*inputs.PatchDiam1
              % Isolated component, do nothing
            else
              % Join to "Trunk"
              % 找到最近点在comp的位置
              I = comp(IC);
              % 找到最近点在trunk的位置 
              J = trunk(IT);
              Other(comp) = false;
              % 属于trunk内
              Trunk(comp) = true;
              NonClassified(i) = false;
              % make the connection
              % 添加彼此，互相定义为邻居
              Nei{I} = [Nei{I}; J];
              Nei{J} = [Nei{J}; I];
            end
          end
        end
      end
    end
    % 若NonClassified内仍有集合时，开始下一圈循环
    first_round = false;
    % If "Other" has decreased, do another check with same "distance"
    % npre是未连结的点个数
    if nnz(NonClassified) < npre
        again = true;
        % 更新 新的未连结的点个数
        npre = nnz(NonClassified);
    else
        again = false;
    end
  end
  % k 是额外向外扩展的搜索半径
  k = k+k0; % increase the cell size of the nearby search space
  % Cmin是若点云内不全是树的点时的内圈半径
  Cmin = 3*Cmin; % increase the acceptable component size
end
Forb(Base) = false;
% 更新 邻居关系
cover.neighbor = Nei;

end % End of function

