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

function cover = cover_sets(P,inputs,RelSize)

% ---------------------------------------------------------------------
% COVER_SETS.M          Creates cover sets (surface patches) and their
%                       neighbor-relation for a point cloud
%
% Version 2.00
% Latest update     16 Aug 2017
%
% Copyright (C) 2013-2017 Pasi Raumonen
% ---------------------------------------------------------------------

% Covers the point cloud with small sets, which are along the surface, 
% such that each point belongs at most one cover set; i.e. the cover is 
% a partition of the point cloud. 
% 
% The cover is generated such that at first the point cloud is covered 
% with balls with radius "BallRad". This first cover is such that 
% 1) the minimum distance between the centers is "PatchDiam", and 
% 2) the maximum distance from any point to nearest center is also "PatchDiam".
% Then the first cover of BallRad-balls is used to define a second cover:
% each BallRad-ball "A" defines corresponding cover set "B" in the second cover
% such that "B" contains those points of "A" that are nearer to the center of
% "A" than any other center of BallRad-balls. The BallRad-balls also define 
% the neighbors for the second cover: Let CA and CB denote cover sets in 
% the second cover, and BA and BB their BallRad-balls. Then CB is 
% a neighbor of CA, and vice versa, if BA and CB intersect or 
% BB and CA intersect.
%
% Inputs: 
% P         Point cloud
% inputs    Input stucture, the following fields are needed:
%   PatchDiam   Minimum distance between centers of cover sets; i.e. the
%                   minimum diameter of a cover set. If "RelSize" is given
%                   as input, then there is minimum and maximum PatchDiam
% 	BallRad     Radius of the balls used to generate the cover sets, these 
%                   balls are also used to determine the neighbors and the 
%                   cover set characteristics              
%   nmin        Minimum number of points in a rcov-ball
% RelSize   Relative cover set size for each point
%
% Outputs:
% cover     Structure array containing the followin fields:
%   ball        Cover sets, (n_sets x 1)-cell
%   center      Center points of the cover sets, (n_sets x 1)-vector
%   neighbor    Neighboring cover sets of each cover set, (n_sets x 1)-cell

if ~isa(P,'double')
    P = double(P);
end

%% Large balls and centers

% np是点云集P中的点数
np = size(P,1);
% new元胞数组(np x 1)
Ball = cell(np,1);  % the large balls, used to generate the cover sets and their neighbors
% new(npx1)全为0的矩阵
Cen = zeros(np,1,'uint32');  % the center points of the balls/cover sets
% new(npx1)全为1的矩阵--布尔类型
NotExa = true(np,1); % the points not yet examined
% new(npx1)全为10^8的矩阵
Dist = 1e8*ones(np,1,'single');  % distance of point to the closest center 
% new(npx1)全为0的矩阵
BoP = zeros(np,1,'uint32');  % the balls/cover sets the points belong
nb = 0;             % number of sets generated

% nargin判断输入变量个数
if nargin == 2
    %% Same size cover sets everywhere
    BallRad = inputs.BallRad1;
    PatchDiamMax = inputs.PatchDiam1;
    nmin = inputs.nmin1;
    % Partition the point cloud into cubes for quick neighbor search
    % partition是一个复杂的元胞数组
    % CC是P内每个点位于三维的正方形位置
    [partition,CC] = cubical_partition(P,BallRad);
    
    % Generate the balls
    Radius = BallRad^2;
    MaxDist = PatchDiamMax^2;
    % randperm生成从1-np不重复随机排列的行向量
    RandPerm = uint32(randperm(np)); % random permutation of points, 
                                     % results in different covers for same input
    for i = 1:np
        if NotExa(RandPerm(i))
            % 任选矩阵内的点Q
            Q = RandPerm(i);
            % 取Q点所在的正方体周围共27个正方体为研究对象
            points = partition(CC(Q,1)-1:CC(Q,1)+1,CC(Q,2)-1:CC(Q,2)+1,CC(Q,3)-1:CC(Q,3)+1);
            % vertcat内置函数, 将points平凑为n行1列的矩阵
            points = vertcat(points{:});
            % 对研究对象进行归一化，以Q点所在正方体为圆心
            V = [P(points,1)-P(Q,1) P(points,2)-P(Q,2) P(points,3)-P(Q,3)];
            % sum(a,2)对矩阵a的行求和,sum(a,1)对矩阵a的列求和
            % 各个点对中心点的距离
            dist = sum(V.*V,2);
            % J矩阵为布尔矩阵, 若距离<半径则为1,否则为0
            J = dist < Radius;
            % nnz=number of nonzero内置函数,返回非零的个数
            % 若以Radius为半径的球体内，有足够多的点则进入if
            if nnz(J) >= nmin
                % I为J中不为零位置，所对应的点云坐标点
                I = points(J);
                % 取出dist内J=1的距离值
                d = dist(J);
                % J为距离小于最大距离的值
                J = (dist < MaxDist);
                % 将本次循环内，所有的近邻点，设置为false=0，剔除掉
                NotExa(points(J)) = false;
                % nb = 生成的sets数
                nb = nb+1;
                % 将本次Q点附近的点保存在Ball中
                Ball{nb} = I;
                % 将本次Q点中心点保存在Cen中
                Cen(nb) = Q;
                % D 为取出 len(I)个 10^8的值
                D = Dist(I);
                % L为布尔类型，理论上均为1，因为d的距离一定远远小于D(10^8), 除非存在误差点（离中心点非常非常远），则这步的目的是去除误差点
                L = d < D;
                % I = I 除非存在误差点，否则可以理解为 变量没有变化
                I = I(L);
                % Dist = d 除非存在误差点，否则可以理解为 变量没有变化
                Dist(I) = d(L);
                % 将BoP(就是原始点云P)中，I位置的这些点，归为同一类(nb的值)，聚类，并给他们类的序号
                BoP(I) = nb;
            end
        end
    end
    % 测试点云1.txt中，共1910133个点，共14853个类别（每次运行结果会有区别，因为是随机取的中心点）
    % 只保留nb个类别的中心周围点，，，（ , 后的 : 没看懂）
    Ball = Ball(1:nb,:);
    % 只保留nb个类别的中心点
    Cen = Cen(1:nb);
else
    %% Use relative sizes (the size varies)
    % Partition the point cloud into cubes
    BallRad = inputs.BallRad2;
    PatchDiamMin = inputs.PatchDiam2Min;
    PatchDiamMax = inputs.PatchDiam2Max;
    nmin = inputs.nmin2;
    MRS = PatchDiamMin/PatchDiamMax;
    r = double(1.5*(double(min(RelSize))/256*(1-MRS)+MRS)*BallRad+1e-5); % minimum radius
    NE = 1+ceil(BallRad/r);
    if NE > 4
        r = PatchDiamMax/4;
        NE = 1+ceil(BallRad/r);
    end
    [Partition,CC,~,Cubes] = cubical_partition(P,r,NE);
    
    I = RelSize == 0; % Don't use points with no size determined
    NotExa(I) = false;
    
    % Define random permutation of points (results in different covers for same input)
    % so that first small sets are generated
    RandPerm = zeros(np,1,'uint32');
    I = RelSize <= 32;
    ind = uint32(1:1:np)';
    I = ind(I);
    t1 = length(I); 
    RandPerm(1:1:t1) = I(randperm(t1));
    I = RelSize <= 128 & RelSize > 32;
    I = ind(I);
    t2 = length(I);
    RandPerm(t1+1:1:t1+t2) = I(randperm(t2));
    t2 = t2+t1;
    I = RelSize > 128;
    I = ind(I);
    t3 = length(I);
    RandPerm(t2+1:1:t2+t3) = I(randperm(t3));
    clearvars ind I
    %RandPerm = uint32(randperm(np)); % random permutation of points

    Point = zeros(round(np/1000),1,'uint32');
    e = BallRad-PatchDiamMax;
    for i = 1:np
        if NotExa(RandPerm(i))
            Q = RandPerm(i);
            rs = double(RelSize(Q))/256*(1-MRS)+MRS; % relative radius
            Dmin = PatchDiamMax*rs;
            Radius = Dmin+sqrt(rs)*e;
            N = ceil(Radius/r); % = number of cells needed to include the ball
            cubes = Cubes(CC(Q,1)-N:CC(Q,1)+N,CC(Q,2)-N:CC(Q,2)+N,CC(Q,3)-N:CC(Q,3)+N);
            I = cubes > 0;
            cubes = cubes(I);
            Par = Partition(cubes);
            S = cellfun('length',Par);
            stop = cumsum(S);
            start = [0; stop]+1;
            for k = 1:length(stop)
                Point(start(k):stop(k)) = Par{k};
            end
            points = Point(1:stop(k));
            V = [P(points,1)-P(Q,1) P(points,2)-P(Q,2) P(points,3)-P(Q,3)];
            dist = sum(V.*V,2);
            J = dist < Radius^2;
            if nnz(J) >= nmin
                I = points(J);
                d = dist(J);
                J = (dist < Dmin^2);
                NotExa(points(J)) = false;
                nb = nb+1;
                Ball{nb} = I;
                Cen(nb) = Q;
                D = Dist(I);
                L = d < D;
                I = I(L);
                Dist(I) = d(L);
                BoP(I) = nb;
            end
        end
    end
    Ball = Ball(1:nb,:);
    Cen = Cen(1:nb);
end
clearvars RandPerm NotExa Dist

%% Cover sets
% Number of points in each ball and index of each point in its ball

% 全为0的, nbx1矩阵, nb为聚合的类数
Num = zeros(nb,1,'uint32');
% 全为0的, npx1矩阵, np为点云总点数
Ind = zeros(np,1,'uint32');
% 下面的循环 用来统计 每个类别的点数
for i = 1:np
    if BoP(i) > 0
        % 该类别的个数 + 1
        Num(BoP(i)) = Num(BoP(i))+1;
        % 为后面Bal的索引位置
        Ind(i) = Num(BoP(i));
    end
end

% Initialization of the "Bal"
% 新建一个元胞数组, nb为类别数
Bal = cell(nb,1);
for i = 1:nb
    % 每个元胞数组对应的元素位, 初始化相对应长度(类别的点数)的矩阵
    Bal{i} = zeros(Num(i),1,'uint32');
end

% Define the "Bal"
for i = 1:np
    if BoP(i) > 0
        % 元胞数组Bal中, 第BoP(i)的类别, 的第Ind(i)个位置的元素, 赋值为i(在点云数组中的排序位置,第几个)
        Bal{BoP(i),1}(Ind(i)) = i;
    end
end

%% Neighbors
% 定义 "邻居"
% A\B两个集合,若A集合中,存在B集合内的点,那么A\B两个集合是邻居
% Define neighbors. Sets A and B are neighbors if the large ball of A 
% contains points of B. Notice that this is not a symmetric relation.
Nei = cell(nb,1);
Fal = false(nb,1);

% ********
% Ball元胞数组内保存的是每次的"近邻点"
% BoP内保存的是最终每个点所属的类别
% 因为存在,某个点同时被划分为多个类别的情况
% 因此,下面的for循环的目的,就是通过对比
% Ball内每次循环的结果,与最终的结果,的区别,从而找到"近邻点",即两次所属类别不同的点
% ********   

for i = 1:nb
    % B 中是,中心点附近的近邻点
    B = Ball{i};        % the points in the big ball of cover set "i" 
    % I 中是,多个中心点附近的近邻点(布尔类型)
    I = (BoP(B) ~= i);
    % N 中是,从 B 中取出 I 所对应的点 在点云集合P内的位置  
    N = B(I);           % the points of B not in the cover set "i"
     % N 中是,从 BOP 取出 I 所对应的点 的 "邻居集合"
    N = BoP(N);
    
    % select the unique elements of N:
    % 返回N的长度/数量
    n = length(N);
    if n > 2
        % Include 是 nx1 全为1布尔类型的矩阵
        Include = true(n,1);

        for j = 1:n
            % 
            if ~Fal(N(j))
                % 如果是邻居集合,那么Fal中为true, 否则是false
                Fal(N(j)) = true;
            else
                % 如果是1,则是一个新的邻居集合开始,相同邻居集合的点都为0,直到出现新的1,代表新的邻居集合开始
                Include(j) = false;
            end
        end
        % 初始化,方便下次筛选
        Fal(N) = false;
        % N 取出 本轮该中心点(该集合)的所有邻居集合的 序号
        N = N(Include);
    elseif n == 2
        if N(1) == N(2)
            N = N(1);
        end
    end
    
    % 保存在Nei元胞数组中
    Nei{i} = uint32(N);
end

% Make the relation symmetric by adding, if needed, A as B's neighbor 
% in the case B is A's neighbor

% 下面循环的作用
% A是B的邻居,那么B也是A的邻居
% 若A中有B,但是B中没A,则把A添加到B中

for i = 1:nb
    N = Nei{i};
    for j = 1:length(N)
        K = (Nei{N(j)} == i);
        % 若K内均为0,则返回0;否则只要均为非0,则返回1
        if ~any(K)
            % 如果 邻居集合中没有自己,那么就把自己添加到邻居集合内
            Nei{N(j)} = uint32([Nei{N(j)}; i]);
        end
    end
end

% Define output
clear cover
% Bal是分好的类别
cover.ball = Bal;
% Cen是每个类别对应的中心点
cover.center = Cen;
% Nei是每个类别自己的邻居类别
cover.neighbor = Nei;

%% Display statistics
%disp(['    ',num2str(nb),' cover sets, points not covered: ',num2str(np-nnz(BoP))])