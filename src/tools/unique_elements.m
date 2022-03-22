function Set = unique_elements(Set,False)
% 独一无二的变量
% 将set内重复的值删掉


n = length(Set);
if n > 2
    I = true(n,1);
    for j = 1:n
        if ~False(Set(j))
            False(Set(j)) = true;
        else
            I(j) = false;
        end
    end
    Set = Set(I);
elseif n == 2
    if Set(1) == Set(2)
        Set = Set(1);
    end
end