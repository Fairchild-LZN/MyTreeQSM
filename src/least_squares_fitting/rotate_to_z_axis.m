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

function [R,D,a] = rotate_to_z_axis(Vec)

% 返回被旋转的矩阵
% 生成一个和圆柱的z轴方向相同的
% 但垂直与单位向量的向量

% --------------------------------------------------------------------------
% ROTATE_TO_Z_AXIS.M   Forms the rotation matrix to rotate the vector to 
%                           a point along the positive z-axis. 
%
% Input 
% Vec      Vector (3 x 1)
%
% Output 
% R        Rotation matrix, with R * Vec = [0 0 z]', z > 0 

% 叉积
% 返回一个同时和A、B向量垂直的向量
D = cross(Vec,[0 0 1]);
% 若D不是单位向量
if norm(D) > 0
    % 求z轴的反余弦
    a = acos(Vec(3));
    R = rotation_matrix(D,a);
else
    R = eye(3);
    a = 0;
    D = [1 0 0];
end
