function B = positive(A)
% 如果是正值，则不变
% 如果是负值，则使其等于0

%% 调用函数
% fun1.m
% main.m

%% 功能部分
    B=max(0,A);
    %{
    if A >= 0
        B = A;
    else
        B = 0;
    end
    %}
end