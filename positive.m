function B = positive(A)
% �������ֵ���򲻱�
% ����Ǹ�ֵ����ʹ�����0

%% ���ú���
% fun1.m
% main.m

%% ���ܲ���
    B=max(0,A);
    %{
    if A >= 0
        B = A;
    else
        B = 0;
    end
    %}
end