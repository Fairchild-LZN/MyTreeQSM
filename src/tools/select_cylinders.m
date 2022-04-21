function cylinder = select_cylinders(cylinder,Ind)

% cylinder元胞数组中的每行名字
Names = fieldnames(cylinder);
% 元胞数组的变量个数
n = size(Names,1);
% 选择cylinder中第Ind位置的所有信息
for i = 1:n
    cylinder.(Names{i}) = cylinder.(Names{i})(Ind,:);
end