clear all;
clc;
% ��

%% *************************����ʽ��*****************************
% W_lasti_i_t = D_hpX_i_t + D_Xmax_i_t + D_lpX_i + sigma_max_Cj + ...
% (num(Pi)-1)*sl - positive(sigma_delta_h_i_t - t) - Ci*(1+max(alfa_neg_X / alfa_pos_X))

% Ri = max(W_lasti_i_t - t + Ci)


%% ************************��չʾ�������˽ṹ��*******************
%map=imread('topology_eg.jpg');
%imshow(map);


%% *****************������ģ�Ͷ��崫�������***********************
%disp('��������ֱ�Ӽ��롣');

%Tao = input('\n�������������������\n�� = ');
Tao = xlsread('data',1,'A2');

%R = input('\n��������:\nR=');   % Mbit/s
R = xlsread('data',1,'B2');
Alfa = R * [0.5 0.5 0.25 0.75];
global alfa_neg_A;
global alfa_pos_A;
global alfa_neg_B;
global alfa_pos_B;
       alfa_pos_A = Alfa(1);
       alfa_neg_A = -Alfa(2);
       alfa_pos_B = Alfa(3);
       alfa_neg_B = -Alfa(4);

global I;
       %I = input('\n����뱻�о�����ĽǱ�\nI=');  %I�Ǳ��о��������ĽǱ꣬Ҳ���������������е�λ�á�
       I = xlsread('data',1,'C2');
global sl;
       %sl = input('\n�����switch latency\nsl = ');
       sl = xlsread('data',1,'D2');
global t;
       %t=input('\n�����I����ʱ��t\n');
       t = xlsread('data',1,'E2');
%disp('*************************************************');
%disp('ע�⣡');
%disp('�������������Ӧ��Ϊ�����������ʽ��');
disp(' ');

disp('***********************������ʱ��C��***********************');
global C;
       for i = 1:Tao
           %disp(['�����ÿһ���������ĵ�λ֡����ʱ��C',num2str(i),'/��',num2str(Tao),'����ms��']);
           %C(i) = input('C = ');
           %eval(['C(i) = xlsread(','''data''',num2str(1),'''F',num2str(i+1),''')']);
           eval(['C(i) = xlsread(''data'',1,','''F',num2str(i+1),''');']);
       end
       C
       %C=[40 40 40 40 40 120];                 %����ʱ�䣨����֡����
       %C = input('\n�����ÿһ���������ĵ�λ֡����ʱ��(ms)��\nC=');
       
disp('***********************����������T��***********************');
global T;
       for i = 1:Tao
           %disp(['�����ÿһ���������ĵ�λ֡��������T',num2str(i),'/��',num2str(Tao),'����ms��']);
           %T(i) = input('T = ');
           eval(['T(i) = xlsread(''data'',1,','''G',num2str(i+1),''');']);
       end
       T
       %T=[2000 2000 2000 160 2000 2000];       %��������
       %T = input('\n�����ÿһ���������ĵ�λ֡��������(ms)��\nT=');
       
disp('***********************�����ȼ���***********************');
global Pr;
       for i = 1:Tao
           %disp(['��Ӧ�1����n���μ���������������ȼ�Pr',num2str(i),'/��',num2str(Tao),'�� (1/2/3)']);
           %Pr_temp = input('Pr = ');
           eval(['Pr_temp = xlsread(''data'',1,','''H',num2str(i+1),''');']);
           while (Pr_temp ~= 1)&&(Pr_temp ~= 2)&&(Pr_temp ~= 3)         %�����������ȼ��Ƿ���1/2/3��Χ��
                 disp(['������',num2str(i),'����Ч���ݣ����ȼ�ֻ��Ϊ1/2/3��']);
                 Pr_temp = input('Pr = ');
           end
           Pr(i) = Pr_temp;
       end
       Pr
       %Pr=[2 2 2 3 1 1];                       %���ȼ�
       %Pr = input('\n��Ӧ�1����n��˳�����ÿһ�������������ȼ�(1/2/3)��\nPr=');
       
%�����Class A/B/C��λ�á�
if ismember(3,Pr)       %�����Class A��������
    %{
    ���´��뱻find��������
    A = 1;
    for i=1:length(Pr)
        if Pr(i) == 3       %Class A���ȼ������3
            j_a(A) = i;     %����Pr�е�Class A
        end
            A = A + 1;
    end
    A = 1;
    for i=1:length(j_a)
        if j_a(i) ~= 0
            temp = j_a(i);
            J_A(A) = temp;      %����j_a�еķ�0Ԫ�أ�J_AΪclassA������������ 
            A = A + 1;
        end
    end
    
    A = 1;
    for i = 1:length(j_a)
        if (j_a(i) ~= 0) && (j_a(i) ~= I)       %����J_A�еķ�I������
            temp = j_a(i);
            j_Reala(A) = temp;
        end
        A = A + 1;
    end
    A = 1;
    for i=1:length(j_Reala)
        if j_Reala(i) ~= 0
            temp = j_Reala(i);
            J_RealA(A) = temp;      %����j_Reala�еķ�0Ԫ�أ�J_RealAΪclassA�ĳ����о����������������  //[4]
            A = A + 1;              %��һ������ŵ���if����Ͳ����𵽳�ȥ0Ԫ�ص�����
        end
    end
    %}
    J_A = find(Pr == 3);
    J_A_temp = find(Pr == 3);
    %if ismember(I,J_A)      %���I��Class A
    J_A_temp(find(J_A_temp == I)) = [];   %ɾ��I���ڵ�λ�õ�����
        %location_I = find(J_A == I);
        %J_A(location_I) = [];
    J_RealA = J_A_temp;      %J_RealA�������
    %end
    clear location_I J_A_temp;
end

if ismember(2,Pr)       %�����Class B��������
    %{
    ���´��뱻find��������
    B = 1;
    for i=1:length(Pr)
        if Pr(i) == 2       %Class B���ȼ���2
            j_b(B) = i;     %����Pr�е�Class B
        end
        B = B + 1;
    end
    B = 1;
    for i=1:length(j_b)
        if j_b(i) ~= 0
            temp = j_b(i);
            J_B(B) = temp;      %����j_b�еķ�0Ԫ�أ�J_BΪclassB������������  //[1 2 3]
            B = B + 1;
        end
    end
    
    B = 1;
    for i = 1:length(j_b)
        if (j_b(i) ~= 0) && (j_b(i) ~= I)
            temp = j_b(i);
            j_Realb(B) = temp;  %ȥ�������������ͬClass�ı��о�����I
        end
        B = B + 1;
    end
    B = 1;
    for i=1:length(j_Realb)
        if j_Realb(i) ~= 0
            temp = j_Realb(i);
            J_RealB(B) = temp;      %����j_Realb�еķ�0Ԫ�أ�J_RealBΪclassB�ĳ����о����������������  //[2 3]
            B = B + 1;              %��һ������ŵ���if����Ͳ����𵽳�ȥ0Ԫ�ص�����
        end
    end
    %}
    J_B = find(Pr == 2);
    J_B_temp = find(Pr == 2);
    %if ismember(I,J_B)      %���I��Class B
    J_B_temp(find(J_B_temp == I)) = [];
        %location_I = find(J_B == I);
        %J_B(location_I) = [];
    J_RealB = J_B_temp;
    %end
    clear location_I J_B_temp;
end

%I���Լ��޷���������
S_first_ii_max_i = 0;
S_last_ii_min_i = 0;
A_ii = 0;


%% ************����S_InŪ����ά���󲢵õ�S_InSeq��*****************
%ÿһ�д���һ������˿ڣ�ÿһ��û��ʵ�ʺ��壨һ���е�ĳһ�д�������˿��е�һ��֡��������ÿһҳ����һ��������

%Num = input('\n������������\nNum = ');
Num = xlsread('data',1,'I2');

%Max_Input = input('\n��༸������˿ڣ�\nMax_Input = ');
Max_Input = xlsread('data',1,'J2');

%Max_InputNum = input('\nĳ������˿ڵģ������������\nMax_InputNum = ');      %��������ÿһ��Ԫ������
Max_InputNum = xlsread('data',1,'K2');
%********************************************************************************
%********************************************************************************
%********************************************************************************
%*****************************������In.xls��*************************************
global S_InSeq;               %������������
S_In = zeros(Max_Input,Max_InputNum,Num);  %��Max_Input������༸������˿ڣ�
                                           %��Max_InputNum����ĳ������˿���������������
                                           %ҳNum������������

judge = exist('In.xls','file'); %�ж��Ƿ���ھɵ�In.xls
if judge == 2   %���Ŀǰ��������In.xls�ļ�
    select = input('�Ƿ���Ҫɾ������In.xls���������ɣ� Y/N\n','s');
    if select == 'Y'||'y'    %��Ҫɾ��
        delete('In.xls');   %ɾ��In.xls�Է������������������In.xls���
        disp('�ѳɹ�ɾ���ɱ�');
        %����In.xls����Լ�������
        disp('���ڴ����±���');
        for i = 1:Num
            eval(['xlswrite(''In'',S_In(:,:,',num2str(i),'),',num2str(i),');']);    %Sheet i�ʹ���Si����������.�Ѵ�0����д��In.xls��
            switch_loca = {'������λ��'};
            eval(['xlswrite(''In'',switch_loca,',num2str(i),',''A',num2str(Max_Input + 2),''');']);    %����������е�λ�ô����������뽻����λ�õĸ���
            eval(['xlswrite(''In'',',num2str(1),',',num2str(i),',''A',num2str(Max_Input + 3),''');']);  %Ԥ��ÿ��λ�ö���1
            
            title1 = {'Sheet'};
            title2 = {'����S'};
            eval(['xlswrite(''In'',title1,',num2str(i),',''A',num2str(Max_Input + 5),''');']);
            eval(['xlswrite(''In'',',num2str(i),',',num2str(i),',''B',num2str(Max_Input + 5),''');']);
            eval(['xlswrite(''In'',title2,',num2str(i),',''C',num2str(Max_Input + 5),''');']);
            eval(['xlswrite(''In'',',num2str(i),',',num2str(i),',''D',num2str(Max_Input + 5),''');']);
        end
    end
else if judge == 0       %���������In.xls
        %����In.xls����Լ�������    
        disp('��������In.xls��񡭡�');
        for i = 1:Num
            eval(['xlswrite(''In'',S_In(:,:,',num2str(i),'),',num2str(i),');']);    %Sheet i�ʹ���Si����������.�Ѵ�0����д��In.xls��
            switch_loca = {'������λ��'};
            eval(['xlswrite(''In'',switch_loca,',num2str(i),',''A',num2str(Max_Input + 2),''');']);    %����������е�λ�ô����������뽻����λ�õĸ���
            eval(['xlswrite(''In'',',num2str(1),',',num2str(i),',''A',num2str(Max_Input + 3),''');']);  %Ԥ��ÿ��λ�ö���1

            title1 = {'Sheet'};
            title2 = {'����S'};
            eval(['xlswrite(''In'',title1,',num2str(i),',''A',num2str(Max_Input + 5),''');']);
            eval(['xlswrite(''In'',',num2str(i),',',num2str(i),',''B',num2str(Max_Input + 5),''');']);
            eval(['xlswrite(''In'',title2,',num2str(i),',''C',num2str(Max_Input + 5),''');']);
            eval(['xlswrite(''In'',',num2str(i),',',num2str(i),',''D',num2str(Max_Input + 5),''');']);
        end
    end
end


disp('���In.xls�ļ�,����Sheet_X�м��뽻����S_X�����ݣ��������رձ�񣬲��û��س���');
disp('�������趯�������س���');
disp(' ');
pause;
disp('----------------------------');
disp('���н������������ݣ�������ɣ�');
disp('���Ժ󡭡����ڴ�����');
disp('----------------------------');
disp(' ');

%��In.xls����ڵ����ݶ���S_In�ڡ�Sheet1ΪS1��Sheet2ΪS2
%�ɰ棺
%{
 ����Ϊ�ɰ����������� 
       for i = 1:Num
           disp('********************');
           for j = 1:Max_Input
               disp(['�����S',num2str(i),'�������ĵ�',num2str(j),'������˿ڵ����ݣ�']);
               disp(['���ݼ����ʽΪ����������֤��Ԫ�ظ�����������������������磺[a b c]��']);
               S_In(j,:,i) = input(' ');
               disp(' ');
           end
       end
%}  
%�°棺
for i = 1:Num
    for j = 1:Max_Input
            range = [num2str(j),':',num2str(j)];    %��ȡ��j��
            S_In(j,:,i) = xlsread('In',i,range);     %��һ��һ�ж�,�õ�����S_In
    end
end


%S1=[1 0 0;2 3 6;0 0 0];
%S2=[1 2 3;4 0 0;5 0 0];
%S3=[2 0 0;3 0 0;6 0 0];

%ʵ�֣���1ҳ��S3����2ҳ��S1����3ҳ��S2���õ�S_InSeq
%�ɰ棺
%{
�ɰ������˳��
for i = 1:Num
    disp(['���ǵ�',num2str(i),'������������ȷ������S��������������ִ���x��S_In(:,:,x)']);
    x=input('');
    temp = S_In(:,:,x);
    S_InSeq(:,:,i) = temp;       %S_InSeq�Ǿ�������˳��Ľ������ṹ�����һҳ�ʹ����һ��������
    disp(' ');
end
%}
%�°棺
for i = 1:Num
    range = [num2str(1),':',num2str(Max_Input)];
    temp(:,:) = xlsread('In',i,range);     %��ȡ��S1��S2��S3������
    seq = xlsread('In',i,['A',num2str(Max_Input + 3)]);     %seq���������������λ�ã��ڼ�����
    S_SwitchSeq(i) = seq;   %Ū��S���ǵڼ�����������
    S_InSeq(:,:,seq) = temp;
end

for i = 1:Num
    n = 1;
    while ~isequal(S_In(:,:,n),S_InSeq(:,:,i))
        n = n + 1;
    end
    S_ReSwitchSeq(i) = n;   %Ū���ڼ�����������S����
end

%********************************************************************************
%********************************************************************************
%********************************************************************************
%********************************************************************************

disp('***********************����������ֲ�����ͼ��***********************');
S_InSeq
clear judge select switch_loca title1 title2 range temp seq i j n;


%% **************������S_Out��S_OutSeq��***********************
%Max_Output = input('\n��༸������˿ڣ�\nMax_Output = ');
Max_Output = xlsread('data',1,'L2');

%Max_OutputNum = input('\nĳ������˿ڵģ�����������\nMax_OutputNum = ');     %��������ÿһ��Ԫ������
Max_OutputNum = xlsread('data',1,'M2');
disp(['���о���Ԫ�ظ���Ϊ',num2str(Max_OutputNum),'��']);

global S_OutSeq;               %�����������
S_Out = zeros(Max_Output,Max_OutputNum,Num);  %�д�����༸������˿ڣ�
                                              %�д���ĳ������˿��������������
                                              %ҳ������������
   for i = 1:Num        %Ϊ�˵õ�ÿһ�����������������
       disp('*******************************************');
       for j = 1:Max_Output
           disp(' ');
           disp(['�����S',num2str(i),'��������������ݣ�']);
           disp('���ݼ����ʽΪ����������֤��Ԫ�ظ��� = ���������������磺[a b c]��');
           S_Out(j,:,i) = input(' ');
           disp(' ');
       end
   end
clear i j;
%ʵ�֣���1ҳ��S3����2ҳ��S1����3ҳ��S2��
for i = 1:Num
    disp(['���ǵ�',num2str(i),'������������ȷ������S��������������ִ���x��S_Out(:,:,x)']);
    x=input('');
    temp = S_Out(:,:,x);
    S_OutSeq(:,:,i) = temp;       %S_OutSeq�Ǿ�������˳��Ľ������ṹ�����һҳ�ʹ����һ��������
    S_OutSeqNo(i) = x;            %S_OutSeqNo����ڼ�����������S��������ֱ�ӵõ����
    disp(' ');
end
clear x temp i;
%�����S_Out����ΪҪ����ֱ����S����ȷ�������������


%% ***************��Ѱ��Pi����I��·����************************
%·������ĵ�һ�д���������ţ�ÿһ�д���һ���������ϵ�����˿�
p = 0;
for i = 1:Num
    for j = 1:Max_Output
        m = 0;
        for k = 1:Max_OutputNum
            if S_OutSeq(j,k,i) == I     %�б��о�������ڣ���ζ���������˿�����Pi
                m=m+1;
                p=p+1;
            end
            if m == 1
                Pi(:,p) = [S_OutSeqNo(i),j];        %PiΪ·�����󣨵�һ���ǽ�������ţ�
            end
        end
    end
end
clear i j m p;


%% ************************��ռλ��1��**********************************
%{



































%}


%% **********************���ҳ�source_i/j��************************************
%�����ҵ�ÿһ���������ĳ�����λ�ã��Խ�����Ϊ��㣬����2������������S3����һ������������¼�ľ���1��
source = zeros(Tao,Num);    %һ����Tao�������������Num��������Ľ�����������Ϊ�˷�ֹ��
                            %����һ��������ֻ������һ������������һ������
count_source = 1;
for i = 1:Num
    for j = 1:Max_Input
        n = 0;
        for k = 1:Max_InputNum
            if S_In(j,k,i) ~= 0
                n = n + 1;
            end
        end
        if n == 1    %��һ������˿���ֻ��һ�����롣
            loca_source = find(S_In(j,:,i) ~= 0);    %�ҵ�ÿһ������˿��еĵ��������λ�ã�������Ϊ����������ԣ�
            source_flow = S_In(j,loca_source,i);     %����˿��еĵ�������ĽǱ�
            count_source = 1;
            while ~isequal(S_In(:,:,i),S_InSeq(:,:,count_source))     
                count_source = count_source + 1;      
            end
            source(source_flow,i) = count_source;     %source�������壺
                                                      %��ֵ�����n������������������������ţ�
                                                      %���������������Ǳ�
        end
    end
end
for i = 1:Tao
    temp = find(source(i,:) ~= 0);
    Flow_Source(i) = source(i,temp);    %Flow_Source������ÿһ������������ʼλ��
end
clear source source_flow loca_source count_source i temp;


%% **********************���ҳ�first_ij��*************************************
%���������ɨ��һ�飬
%�Ѿ���i��j�Ľ�������S_InSeq���Ÿ��򣬰ѵ�һ���ó�������first_ij�����һ������last_ij
disp(' ');
disp('****************************');
disp('�����¼����first_ij��last_ij');
%�ҵ�����Class B��������
if ismember(2,Pr)
    judge_first = 0;                %�ж��Ƿ�Ϊfirst_ij���жϲ���
    for J = 1:length(J_B)
        for i = 1:length(Pi)        %��Pi�ĵ�һ�е�Pi�����һ��
            %��ҪŪ��i=1ʱ���ĸ���������
            Switch_No = Pi(1,i);    %Pi�ϵĵ�i����������S��Switch_No����Ҳ����S_Out�еĵڡ�Switch_No��ҳ
            for j = 1:Max_Output
                m = 0;
                %{
                ���������������Ϊ�޷��ų���first_11��I = J_B(J)��������ء�
                if (ismember(I,S_Out(j,:,i))) && (ismember(J_B(J),S_Out(j,:,i)))
                    m=m+1;
                end
                %}
                if ismember(J_B(J),S_Out(j,:,Switch_No))
                    for k = 1:Max_OutputNum
                        if (S_Out(j,k,Switch_No) == I) || (S_Out(j,k,Switch_No) == J_B(J))
                            m=m+1;
                        end
                    end
                end
                if m == 2                % m==2˵���������м���I��Ҳ��J_B(J),��ѭ���ֵ��ĸ�����
                    %cmpt(n) = i;         %(cmpt��compete����˼������ڼ���������)
                    %cmpt_out(n) = j;     %(cmpt_out����ĳ���������ϣ��������ڵĵڼ�������˿�)
                    %n = n + 1;
                    m = 0;
                    judge_first = judge_first + 1;      %judge_first����1
                end
            end
            if judge_first == 1;                        %�ڵ�һ�����ӵ�ʱ��Ҳ���ǵ�һ������м���IҲ��J_B��ʱ����first_ij
                first_ij = Switch_No;
                for l = 1:Max_Output
                    if ismember(I,S_Out(l,:,first_ij))
                        first_ij_port = l;
                    end
                end
                disp(['first_',num2str(I),num2str(J_B(J)),' = S',num2str(first_ij),num2str(first_ij_port),'��������S',num2str(first_ij)]);
                %���ýṹ��Cmpt.first���洢first_ij
                %��ѯfirst_ij
                %������I = J_B(J)�����
                Cmpt.first(I,J_B(J),1) = first_ij;  %�д���I���д���j����һҳ����������ţ��ڶ�ҳ����˿ں�
                Cmpt.first(I,J_B(J),2) = first_ij_port;
                break
            end
        end
        judge_first = 0;        %��Pi����һ��֮��judge_first���㣬����һ�����ȼ��ĸ���
    end
    clear i j k m first_ij_port l;
end

%�����е�Class A��
if ismember(3,Pr)
    judge_first = 0;                %�ж��Ƿ�Ϊfirst_ij���жϲ���
    for J = 1:length(J_A)
        [hang lie] = size(Pi);
        for i = 1:lie        %��Pi�ĵ�һ�е�Pi�����һ��
            %��ҪŪ��i=1ʱ���ĸ���������
            Switch_No = Pi(1,i);    %Pi�ϵĵ�i����������S��Switch_No����Ҳ����S_Out�еĵڡ�Switch_No��ҳ
            for j = 1:Max_Output
                m = 0;
                if ismember(J_A(J),S_Out(j,:,Switch_No))
                    for k = 1:Max_OutputNum
                        if (S_Out(j,k,Switch_No) == I) || (S_Out(j,k,Switch_No) == J_A(J))
                            m=m+1;
                        end
                    end
                end
                if m == 2                % m==2˵���������м���I��Ҳ��J_A(J),��ѭ���ֵ��ĸ�����
                    %cmpt(n) = i;         %(cmpt��compete����˼������ڼ���������)
                    %cmpt_out(n) = j;     %(cmpt_out����ĳ���������ϣ��������ڵĵڼ�������˿�)
                    %n = n + 1;
                    m = 0;
                    judge_first = judge_first + 1;
                end
            end
            if judge_first == 1;
                first_ij = Switch_No;
                for l = 1:Max_Output
                    if ismember(I,S_Out(l,:,first_ij))
                        first_ij_port = l;
                    end
                end
                disp(['first_',num2str(I),num2str(J_A(J)),' = S',num2str(first_ij),num2str(first_ij_port),'��������S',num2str(first_ij)]);
                %���ýṹ��Cmpt.first���洢first_ij
                Cmpt.first(I,J_A(J),1) = first_ij;  %�д���I���д���j����һҳ����������ţ��ڶ�ҳ����˿ں�
                Cmpt.first(I,J_A(J),2) = first_ij_port;
                break
            end
        end
        judge_first = 0;
    end
    clear i j k m J judge_first first_ij Switch_No l first_ij_port hang lie;
end


%% *******************���ҳ�last_ij��**************************
%{
������ͨ����ĳ����ͬʱ����I&J���ɴ˽������״̬��
�����߼���ͬʱ���ڣ����ж�Ϊ�棬ͨ����
��������ʧ������һ��������Ϊlast��
Ȼ���ٲ�һ��������һ�����������ĸ��˿ڡ�
���Ǽ���״̬�Ľ�����˳�����ʵ�֣�

��˼·��ֱ������Ѱ��first�ͺ��ˡ�
%}
reverse_Pi = fliplr(Pi);
%�ҵ�����Class B��������
if ismember(2,Pr)
    judge_last = 0;                         %�ж��Ƿ�Ϊlast_ij���жϲ���
    for J = 1:length(J_B)
        for i = 1:length(reverse_Pi)
            Switch_No = reverse_Pi(1,i);    %Pi�ϵĵ�i����������S��Switch_No����Ҳ����S_Out�еĵڡ�Switch_No��ҳ
            for j = 1:Max_Output    %�˿�����
               m = 0;
                if ismember(J_B(J),S_Out(j,:,Switch_No))
                    for k = 1:Max_OutputNum
                        if (S_Out(j,k,Switch_No) == I) || (S_Out(j,k,Switch_No) == J_B(J))
                            m=m+1;
                        end
                    end
                end
                if m == 2     % m==2˵���������м���I��Ҳ��J_B(J),��ѭ���ֵ��ĸ�����
                    %{
                    %����ʼ����last_ij�ļ���״̬������������ͬʱ�������ж�Ϊ�桪��ͨ������
                    %������ͬʱ���ڡ�������һ��������Ϊlast���ٲ�һ������һ�����������ĸ��˿ڡ�
                    whatchdog = 1;      %��������
                    %}
                    m = 0;
                    judge_last = judge_last + 1;
                end
            end
            if judge_last == 1;
                last_ij = Switch_No;
                for l = 1:Max_Output
                    if ismember(I,S_Out(l,:,last_ij))
                        last_ij_port = l;
                    end
                end
                disp(['last_',num2str(I),num2str(J_B(J)),' = S',num2str(last_ij),num2str(last_ij_port),'��������S',num2str(last_ij)]);
                %���ýṹ��Cmpt.last���洢last_ij
                %��ѯlast_ij
                %������I = J_B(J)�����
                Cmpt.last(I,J_B(J),1) = last_ij;
                Cmpt.last(I,J_B(J),2) = last_ij_port;
                break
            end
        end
        judge_last = 0;
    end
    clear J i j m k l judge_last last_ij_port;
end

%�ҵ�����Class A��������
if ismember(3,Pr)
    judge_last = 0;                         %�ж��Ƿ�Ϊlast_ij���жϲ���
    for J = 1:length(J_A)
        [hang lie] = size(reverse_Pi);
        for i = 1:lie
            Switch_No = reverse_Pi(1,i);    %Pi�ϵĵ�i����������S��Switch_No����Ҳ����S_Out�еĵڡ�Switch_No��ҳ
            for j = 1:Max_Output    %�˿�����
                m = 0;
                if ismember(J_A(J),S_Out(j,:,Switch_No))
                    for k = 1:Max_OutputNum
                        if (S_Out(j,k,Switch_No) == I) || (S_Out(j,k,Switch_No) == J_A(J))
                            m=m+1;
                        end
                    end
                end
                if m == 2     % m==2˵���������м���I��Ҳ��J_A(J),��ѭ���ֵ��ĸ�����
                    %{
                    %����ʼ����last_ij�ļ���״̬������������ͬʱ�������ж�Ϊ�桪��ͨ������
                    %������ͬʱ���ڡ�������һ��������Ϊlast���ٲ�һ������һ�����������ĸ��˿ڡ�
                    whatchdog = 1;      %��������
                    %}
                    m = 0;
                    judge_last = judge_last + 1;
                end
            end
            if judge_last == 1;
                last_ij = Switch_No;
                for l = 1:Max_Output
                    if ismember(I,S_Out(l,:,last_ij))
                        last_ij_port = l;
                    end
                end
                disp(['last_',num2str(I),num2str(J_A(J)),' = S',num2str(last_ij),num2str(last_ij_port),'��������S',num2str(last_ij)]);
                %���ýṹ��Cmpt.last���洢last_ij
                %��ѯlast_ij
                %������I = J_A(J)�����
                Cmpt.last(I,J_A(J),1) = last_ij;
                Cmpt.last(I,J_A(J),2) = last_ij_port;
                break
            end
        end
        judge_last = 0;
    end
    clear judge_last J i Switch_No j m k last_ij last_ij_port hang lie;
end
disp('****************************');
disp(' ');


%% ************************��ռλ��2��**********************************
%{















%}


%% ************************��֪����û����***********************
%{
��֪����û����
global S31;
global S11;
global S21;

global SeriaTerm;
global S_S21_min4;
global A_1_4;
global D_classB_1_j_t;
global Rep_1_B;
global delta_S21_1

t=0;
C
%}


%% ********************* ���㡾S_first_ij_min_j���͡�S_first_ij_min_i����ֵ **************************
%{
%����ʾfirst_ij�Ƕ��٣���S11
%disp(['first_ij = S',num2str(CMPT(1,1)),num2str(CMPT(2,1))]);       %CMPT(1,1)��������S����CMPT(2,1)��������˿ڼ�
%�ٸ��first_ij�Ľ������ǵ�No_first_ij������ڡ�2������Ҫ�ҵ�CMPT(1,1)��λ��
%n = 1;
%while ~isequal(S_In(:,:,CMPT(1,1)),S_InSeq(:,:,n))      %��S��CMPT(1,1)�������˽ṹ�еڼ������������
%        n = n + 1;      %���ҵ���first_ij�ǵ�n����������ʱ������ѭ����
%end
%No_first_ij = n;
%}

%S.max = zeros(Max_Output,Tao,Num);   %�������е�S������ҳ����I/J,������ = �ӣ��С������������ƣ�
%S.min = zeros(Max_Output,Tao,Num);   %�С�������˿�������
%����ClassB�����J�ӵ�No_J����������ȥ�ģ���ڡ�1����������J�����ҷ�JԪ�ء�
%***********************************************************
disp('�������J����Class B�ļ��㣺');
count_j = 1;
for J=1:length(J_RealB)
    for i=1:Num                             %ҳ
        for j=1:Max_Input                   %��
            n = 0;
            N = 0;
            for k=1:Max_InputNum            %��
                if S_InSeq(j,k,i) == J_RealB(J)      %һ��һ����J_B
                    n=n+1;
                end
                if S_InSeq(j,k,i) == I      %һ��һ����I
                    N=N+1;
                end
            end
            %********************
            if n == 1                       %������һ�У��������˿ڣ���J����
                m = 0;
                for l = 1:Max_InputNum
                    if (S_InSeq(j,l,i) ~= J_RealB(J)) && (S_InSeq(j,l,i) ~= 0)     %����һ�У��������˿ڣ��з�J����������
                        m=m+1;
                    end
                end
                if m == 0                   %�����ڷ�JԪ�أ�ֻ��J���ڣ�
                    No_J_B(count_j) = i;      %No_J��������ʾS_InSeq�еĵڼ���������
                    count_j = count_j + 1;
                end
            end
            %********************
            if N == 1                       %������һ�У��������˿ڣ���I����
                M = 0;
                for L=1:Max_InputNum
                    if (S_InSeq(j,L,i) ~= I) && (S_InSeq(j,L,i) ~= 0)     %����һ�У��������˿ڣ��з�I����������
                        M=M+1;
                    end

                end
                if M == 0                %�����ڷ�IԪ�أ�ֻ��I���ڣ�
                    No_I = i;
                end
            end
            %********************
        end
    end
end

%����S_first_ij_min_j
%S_first_ij_min_j = C��I��*��No_first_ij - No_J_B + 1��
disp('**************������S_first_ij_min_j��*******************');
disp('********************************************************');
for count_j = 1:length(J_RealB)
    %Name_first_ij = input(['�밴˳�����֮ǰ��¼��first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'�Ľ�������źͶ˿ڣ�������������ʽ���룺\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealB(count_j),1)),S_InSeq(:,:,n))      %��S��CMPT(1,1)�������˽ṹ�еڼ������������
        n = n + 1;      %���ҵ���first_ij�ǵ�n����������ʱ������ѭ����
    end
    %�洢first_ij��λ�á���ţ����֣����˿�
    %��S.min�Ĳ���
    No_first_B(I,J_RealB(count_j)) = n;   
    Name_first_B(I,J_RealB(count_j)) = Cmpt.first(I,J_RealB(count_j),1);
    Port_first_B(I,J_RealB(count_j)) = Cmpt.first(I,J_RealB(count_j),2);
    %No_first_ij = input(['first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'���ڵ�S',...
    %    num2str(Cmpt.first(I,J_RealB(count_j),1)),'�ǵڼ�����������\nNo_first_ij = ']);
    %disp('�����¼����');
    %disp('--------------------------');
    %eval(['S_first_',num2str(I),num2str(J_RealB(count_j)),'_min_',num2str(J_RealB(count_j)),' = ',...
    %    num2str(C(J_RealB(count_j)) * (No_first(I,J_RealB(count_j)) - No_J(count_j) + 1))]);
    disp('--------------------------');
    disp(' ');
    %S.min(Name_first_ij(1),Name_first_ij(2),J_RealB(count_j)) =...
    %    C(J_RealB(count_j)) * (No_first_ij - No_J(count_j) + 1);
    S.min( Cmpt.first(I,J_RealB(count_j),1) , Cmpt.first(I,J_RealB(count_j),2) , J_RealB(count_j)) =...
        C(J_RealB(count_j)) * (No_first_B(I,J_RealB(count_j)) - No_J_B(count_j) + 1);
end
disp('*********************************');
disp(' ');
disp(' ');

%����S_first_ij_min_i
%S_first_ij_min_i = C(I) * (No_first_ij - No_I + 1)
disp('**************������S_first_ij_min_i��*******************');
disp('********************************************************');
disp(['����I�ļ��㣬I = ',num2str(I)]);
if Pr(I) == 3
    disp('I��Class A��������');
else if Pr(I) == 2
        disp('I��Class B��������');
    else
        disp('I��Class C��������');
    end
end
disp(' ');
for count_j = 1:length(J_RealB)
    %Name_first_ij = input(['�밴˳�����֮ǰ��¼��first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'�Ľ�������źͶ˿ڣ�������������ʽ���룺\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealB(count_j),1)),S_InSeq(:,:,n))      
        n = n + 1;      %���ҵ���first_ij�ǵ�n����������ʱ������ѭ����
    end
    No_first_B(I,J_RealB(count_j)) = n;
    %No_first_ij = input(['first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'���ڵ�S',...
    %    num2str(Cmpt.first(I,J_RealB(count_j),1)),'�ǵڼ�����������\nNo_first_ij = ']);
    %disp('�����¼����');
    %disp('--------------------------');
    %eval(['S_first_',num2str(I),num2str(J_RealB(count_j)),'_min_',num2str(I),' = ',...
    %    num2str(C(J_RealB(count_j)) * (No_first(I,J_RealB(count_j)) - No_I + 1))]);
    disp('--------------------------');
    disp(' ');
    %S.min(Name_first_ij(1),Name_first_ij(2),I) =...
    %    C(J_RealB(count_j)) * (No_first_ij - No_I + 1);
    %S.min( Cmpt.first(I,J_RealB(count_j),1) , Cmpt.first(I,J_RealB(count_j),2) , I) =...
    %    C(J_RealB(count_j)) * (No_first_ij - No_I + 1);     %S.min����S_mnΪ��λ�Ĵ洢�ռ䣬��������first/lastΪ��λ�Ĵ洢�ռ�
    S.min( Cmpt.first(I,J_RealB(count_j),1) , Cmpt.first(I,J_RealB(count_j),2) , I) =...
        C(J_RealB(count_j)) * (No_first_B(I,J_RealB(count_j)) - No_I + 1);     %S.min����S_mnΪ��λ�Ĵ洢�ռ䣬��������first/lastΪ��λ�Ĵ洢�ռ�
end
disp('*********************************');
disp(' ');
disp(' ');

%***********************************************************
disp('�������J����Class A�ļ��㣺');
count_j = 1;
for J=1:length(J_RealA)
    for i=1:Num                             %ҳ
        for j=1:Max_Input                   %��
            n = 0;
            N = 0;
            for k=1:Max_InputNum            %��
                if S_InSeq(j,k,i) == J_RealA(J)      %һ��һ����J_A
                    n=n+1;
                end
                if S_InSeq(j,k,i) == I      %һ��һ����I
                    N=N+1;
                end
            end
            %********************
            if n == 1                       %������һ�У��������˿ڣ���J����
                m = 0;
                for l = 1:Max_InputNum
                    if (S_InSeq(j,l,i) ~= J_RealA(J)) && (S_InSeq(j,l,i) ~= 0)     %����һ�У��������˿ڣ��з�J����������
                        m=m+1;
                    end
                end
                if m == 0                   %�����ڷ�JԪ�أ�ֻ��J���ڣ�
                    No_J_A(count_j) = i;      %No_J��������ʾS_InSeq�еĵڼ���������
                    count_j = count_j + 1;
                end
            end
            %********************
            if N == 1                       %������һ�У��������˿ڣ���I����
                M = 0;
                for L=1:Max_InputNum
                    if (S_InSeq(j,L,i) ~= I) && (S_InSeq(j,L,i) ~= 0)     %����һ�У��������˿ڣ��з�I����������
                        M=M+1;
                    end

                end
                if M == 0                %�����ڷ�IԪ�أ�ֻ��I���ڣ�
                    No_I = i;
                end
            end
            %********************
        end
    end
end

%����S_first_ij_min_j
%S_first_ij_min_j = C��I��*��No_first_ij - No_J_A + 1��
disp('**************������S_first_ij_min_j��*******************');
disp('********************************************************');
for count_j = 1:length(J_RealA)
    %Name_first_ij = input(['�밴˳�����֮ǰ��¼��first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'�Ľ�������źͶ˿ڣ�������������ʽ���룺\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealA(count_j),1)),S_InSeq(:,:,n))      %�Ա�S_In��S_InSeq
        n = n + 1;      %���ҵ���first_ij�ǵ�n����������ʱ������ѭ����
    end
    %�洢first_ij��λ�á���ţ����֣����˿�
    %��S.min�Ĳ���
    No_first_A(I,J_RealA(count_j)) = n;   
    Name_first_A(I,J_RealA(count_j)) = Cmpt.first(I,J_RealA(count_j),1);
    Port_first_A(I,J_RealA(count_j)) = Cmpt.first(I,J_RealA(count_j),2);
    disp('--------------------------');
    disp(' ');
    S.min( Cmpt.first(I,J_RealA(count_j),1) , Cmpt.first(I,J_RealA(count_j),2) , J_RealA(count_j)) =...
        C(J_RealA(count_j)) * (No_first_A(I,J_RealA(count_j)) - No_J_A(count_j) + 1);
end
disp('*********************************');
disp(' ');
disp(' ');

%����S_first_ij_min_i
%S_first_ij_min_i = C(I) * (No_first_ij - No_I + 1)
disp('**************������S_first_ij_min_i��*******************');
disp('********************************************************');
disp(['����I�ļ��㣬I = ',num2str(I)]);
if Pr(I) == 3
    disp('I��Class A��������');
else if Pr(I) == 2
        disp('I��Class B��������');
    else
        disp('I��Class C��������');
    end
end
disp(' ');
for count_j = 1:length(J_RealA)
    %Name_first_ij = input(['�밴˳�����֮ǰ��¼��first_',num2str(I),...
    %    num2str(J_RealA(count_j)),'�Ľ�������źͶ˿ڣ�������������ʽ���룺\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealA(count_j),1)),S_InSeq(:,:,n))     
        n = n + 1;      %���ҵ���first_ij�ǵ�n����������ʱ������ѭ����
    end
    No_first_A(I,J_RealA(count_j)) = n;
    disp('--------------------------');
    disp(' ');
    S.min( Cmpt.first(I,J_RealA(count_j),1) , Cmpt.first(I,J_RealA(count_j),2) , I) =...
        C(J_RealA(count_j)) * (No_first_A(I,J_RealA(count_j)) - No_I + 1);     %S.min����S_mnΪ��λ�Ĵ洢�ռ䣬��������first/lastΪ��λ�Ĵ洢�ռ�
end
disp('*********************************');
disp(' ');
disp(' ');
clear i j k l L m M n N count_j J;


%% ***********************������S_first_ij_max_j��S_first_ij_max_i��*******************
S.max(1,1,1) = 40;
S.max(1,1,2) = 360;
S.max(1,1,3) = 40;
S.max(2,1,4) = 40;
S.max(3,1,3) = 40;%��һ����
%{
%% ***********************������S_first_ij_max_j��S_first_ij_max_i��*******************


%S_S11_min_2
S_first_ij_min_j = C(J_b) * (No_first_ij - No_J + 1);
%S_S11_min_1
S_first_ij_min_i = C(I) * (No_first_ij - No_I + 1);

%***********************************************
%����max�д����������TA��
%S_S11_max_2
S_first_ij_max_j = 360;

%S_S11_max_1
S_first_ij_max_i = 40;
%***********************************************

%M_S11_1 = 40
%I��firsti��j�����翪ʼʱ��
M_first_ij_i = C(I) * (No_first_ij - No_I + 1);


%M_ECU1_1 = 0
M_first_i_i = 0;


%M_S21_1 = 80
%��ʾһ��last_ij�Ƕ���
disp(['last_ij = S',num2str(CMPT(1,CMPT1(2))),num2str(CMPT(2,CMPT1(2)))]);
%�ٸ��last_ij�Ľ������ǵ�No_last_ij������ڡ�3������Ҫ�ҵ���CMPT(1,CMPT1(2))����λ��
n = 1;
while ~isequal(S_In(:,:,CMPT(1,CMPT1(2))),S_InSeq(:,:,n))      %�Ƚ����������Ƿ����
    n = n + 1;         %���ҵ���last_ij�ǵ�n����������ʱ������ѭ����
end
No_last_ij = n;        %last_ij�ǵڡ�No_last_ij����������
%����M_last_i_i
M_last_i_i = C(I) * (No_last_ij - No_I+ 1);

%% �� D_hpX_i_t
A_i_j = S_first_ij_min_j - M_first_ij_i;
if Pr(I) == 1
    X = 'Class A';
else if Pr(I) == 2
        X = 'Class B';
    else
        X = 'Class C';
    end
end
disp(['��Ϊ���о�����Ǳ�i=',num2str(I),'�����Ԧ�',num2str(I),'����',X]);
D_classX_IJt = positive(1 + floor((t + S_first_ij_max_i - S_first_ij_min_j + A_i_j) / T(J_b))) * C(J_b);
%}


%% ***********************��M_h_i��*************************
%�����ˣ�����MӦ�ð���first_ij���㡿
%{
%����MӦ����Ҫ����Pi����
for count_m = 1:(length(Pi))     %�������ECU�����ݽ�������ƴ��
    %m_h_i�ǲ�����ECU�����翪ʼʱ��
    m_h_i(I,count_m) = (S_SwitchSeq(Pi(1,count_m)) - Flow_Source(I) + 1) * C(I);     %��һҳ��ECU���ڶ�ҳ��ʼ�ǵ�һ�����������Դ�����
end
%M_h_i�Ǽ���ECU�����翪ʼʱ��
%�����⼸����Ϊ��Ӧ��I����"ĳ����1����"�Ӷ����µ�m_h_i����ά�ȱ仯����
temp = size(m_h_i);
temp1 = zeros(temp(1),1);
temp2 = [temp1,m_h_i]; 
M_h_i = temp2(temp(1),:);
clear count_m m_h_i temp temp1 temp2;
%}

%������B������M
for count_m = 1:length(J_RealB)
    %m_h_i�ǲ�����ECU�����翪ʼʱ��
    %�洢������M_first_ij_i����ʽ������������M_S_mn_i����ʽ
    M_h_i_B(Cmpt.first(I,J_RealB(count_m),1) , Cmpt.first(I,J_RealB(count_m),2) , I) =...
        (S_SwitchSeq(Cmpt.first(I,J_RealB(count_m),1)) - Flow_Source(I) + 1) * C(I);
end
%������A������M
for count_m = 1:length(J_A)
    %m_h_i�ǲ�����ECU�����翪ʼʱ��
    %�洢������M_first_ij_i����ʽ������������M_S_mn_i����ʽ
    M_h_i_A(Cmpt.first(I,J_RealA(count_m),1) , Cmpt.first(I,J_RealA(count_m),2) , I) =...
        (S_SwitchSeq(Cmpt.first(I,J_RealA(count_m),1)) - Flow_Source(I) + 1) * C(I);
end
M_ECU_i = 0;

clear count_m 


%% **************************��A_ij��**************************
%����Class B��A_ij
for i = 1:length(J_RealB)
    % A_ij = S_first_ij_max_j - M_first_ij_i
    % Name_first(I,J_RealB(i))�����˽�������ţ����ƣ���Ҳ����S_mn�е�m
    % Port_first(I,J_RealB(i))����������˿ڣ�Ҳ����S_mn�е�n
    A(I,J_RealB(i)) = S.max(Name_first(I,J_RealB(i)) , Port_first(I,J_RealB(i)) , J_RealB(i)) -...
        M_h_i_B(Name_first(I,J_RealB(i)) , Port_first(I,J_RealB(i)) , I);
end
%����Class A��A_ij
for i = 1:length(J_RealA)
    A(I,J_RealA(i)) = S.max(Name_first_A(I,J_RealA(i)) , Port_first_A(I,J_RealA(i)) , J_RealA(i)) - ...
        M_h_i_A(Name_first_A(I,J_RealA(i)) , Port_first_A(I,J_RealA(i)) , I);
end


%% ***********************��D_ClassX��************************
%ֱ�Ӽ���I���ڵ���һ��
%���Լ���һ���ų�
if ismember(I,J_B)
    for i = 1:length(J_RealB)
        %���I��Class B����I���������ų���
        D_ClassB(I,J_RealB(i)) = positive(1 + floor((t + S.max(Name_first_B(I,J_RealB(i)) , Port_first_B(I,J_RealB(i)) , I) -...
            S.min(Name_first_B(I,J_RealB(i)) , Port_first_B(I,J_RealB(i)) , J_RealB(i)) + A(I,J_RealB(i)))...
            / T(J_RealB(i)))) * C(J_RealB(i));
    end
    %��������I�Լ�
    D_ClassB_I = positive(1 + floor((t + S_first_ii_max_i - S_last_ii_min_i + A_ii)...
        / T(I))) * C(I);
    D_ClassB(I,I) = D_ClassB_I;
    Sum_D_ClassB = sum(sum(D_ClassB));
    clear D_ClassB_I i;
else if ismember(I,J_A)
        for i = 1:length(J_RealA)
            D_ClassA(I,J_RealA(i)) = positive(1 + floor((t + S.max(Name_first_A(I,J_RealB(i)) , Port_first_A(I,J_RealA(i)) , I) -...
                S.min(Name_first_A(I,J_RealA(i)) , Port_first_A(I,J_RealA(i)) , J_RealA(i)) + A(I,J_RealA(i)))...
                / T(J_RealA(i)))) * C(J_RealA(i));
        end
        %��������I�Լ�
        D_ClassA_I = positive(1 + floor((t + S_first_ii_max_i - S_last_ii_min_i + A_ii)...
            / T(I))) * C(I);
        D_ClassA(I,I) = D_ClassA_I;
        Sum_D_ClassA = sum(sum(D_ClassA));
        clear D_ClassA_I i;
    end
end


%% *****************************��Rep_i_X��*******************************
if ismember(I,J_B)
    Rep_B(I) = Sum_D_ClassB * (alfa_neg_B / alfa_pos_B);
else if ismember(I,J_A)
        Rep_A(I) = Sum_D_ClassA * (alfa_neg_A / alfa_pos_A);
    end
end


%% 








%% 

%S_S11_min_2
S_first_ij_min_j = C(J) * (No_first_ij - No_J + 1);
%S_S11_min_1
S_first_ij_min_i = C(I) * (No_first_ij - No_I + 1);


%********************* ���㡾S_first_ij_min_j���͡�S_first_ij_min_i����ֵ **************************
%S_S21_min_2
%S_last_ij_min_j

%S_S21_max_2
%S_last_ij_max_j


%% ***********************************************
%����max�д����������TA��
%S_S11_max_2
S_first_ij_max_j

%S_S11_max_1
S_first_ij_max_i
%***********************************************

%M_S11_1 = 40
%I��firsti��j�����翪ʼʱ��
M_first_ij_i = C(I) * (No_first_ij - No_I + 1);


%M_ECU1_1 = 0
M_first_i_i = 0;


%M_S21_1 = 80
%��ʾһ��last_ij�Ƕ���
disp(['last_ij = S',num2str(CMPT(1,CMPT1(2))),num2str(CMPT(2,CMPT1(2)))]);
%�ٸ��last_ij�Ľ������ǵ�No_last_ij������ڡ�3������Ҫ�ҵ���CMPT(1,CMPT1(2))����λ��
n = 1;
while ~isequal(S_In(:,:,CMPT(1,CMPT1(2))),S_InSeq(:,:,n))      %�Ƚ����������Ƿ����
    n = n + 1;         %���ҵ���last_ij�ǵ�n����������ʱ������ѭ����
end
No_last_ij = n;        %last_ij�ǵڡ�No_last_ij����������
%����M_last_i_i
M_last_i_i = C(I) * (No_last_ij - No_I+ 1);

%% ��D_classX_IJt
A_i_j = S_first_ij_min_j - M_first_ij_i;
if Pr(I) == 1
    X = 'A';
else if Pr(I) == 2
        X = 'B';
    else
        X = 'C';
    end
end
disp(['��Ϊ���о�����Ǳ�i=',num2str(I),'�����Ԧ�',num2str(I),'����Class',X]);
D_classX_IJt = positive(1 + floor((t + S_first_ij_max_i - S_first_ij_min_j + A_i_j) / T(J))) * C(J);

%% ��Rep_i_X


%% �� D_Xmax_i_t


%% �� D_lpX_i


%% �� max(Cj)


%% �� (num(Pi)-1)*sl


%% �� delta_h_i_t


%% �� max(alfa_neg_X / alfa_pos_X)

