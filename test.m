clear all;
clc;

%% *************************【总式】*****************************
%{
  W_lasti_i_t = D_hpX_i_t + D_Xmax_i_t + D_lpX_i + sigma_max_Cj + ...
    (num(Pi)-1)*sl - positive(sigma_delta_h_i_t - t) - Ci*(1+max(alfa_neg_X / alfa_pos_X))

  Ri = max(W_lasti_i_t - t + Ci)
%}


%% ************************【展示网络拓扑结构】*******************
%map=imread('topology_eg.jpg');
%imshow(map);


%% *****************【根据模型定义传输参数】***********************
%{
disp('注意！');
disp('下面输入的数据应该为矩阵或向量形式！');
%Tao = input('\n请键入数据流的总数：\nτ = ');
%}
Tao = xlsread('data',1,'A2');
global I;
I = xlsread('data',1,'C2');  %I是被研究数据流的角标，也就是在数据流其中的位置。
%disp(' ');
%disp('***********************************');
global C;
for i = 1:Tao
           %disp(['请键入每一组数据流的单位帧传输时间C',num2str(i),'/共',num2str(Tao),'个（ms）']);
           %C(i) = input('C = ');
           %eval(['C(i) = xlsread(','''data''',num2str(1),'''F',num2str(i+1),''')']);
           eval(['C(i) = xlsread(''data'',1,','''F',num2str(i+1),''');']);
       end
       C
%{
       for i = 1:Tao
           disp(['请键入每一组数据流的单位帧传输时间C',num2str(i),'/共',num2str(Tao),'个（ms）']);
           C(i) = input('C = ');
           disp(' ');
       end

       %C=[40 40 40 40 40 120];                 %传输时间（代表帧长）
       %C = input('\n请键入每一组数据流的单位帧传输时间(ms)：\nC=');
disp('***********************************');
%}

global T;
for i = 1:Tao
           %disp(['请键入每一组数据流的单位帧发送周期T',num2str(i),'/共',num2str(Tao),'个（ms）']);
           %T(i) = input('T = ');
           eval(['T(i) = xlsread(''data'',1,','''G',num2str(i+1),''');']);
       end
       T
%{
       for i = 1:Tao
           disp(['请键入每一组数据流的单位帧发送周期T',num2str(i),'/共',num2str(Tao),'个（ms）']);
           T(i) = input('T = ');
           disp(' ');
       end
       %T=[2000 2000 2000 160 2000 2000];       %发送周期
       %T = input('\n请键入每一组数据流的单位帧发送周期(ms)：\nT=');
disp('***********************************');
%}
global Pr;
for i = 1:Tao
           %disp(['请从τ1到τn依次键入各数据流的优先级Pr',num2str(i),'/共',num2str(Tao),'个 (1/2/3)']);
           %Pr_temp = input('Pr = ');
           eval(['Pr_temp = xlsread(''data'',1,','''H',num2str(i+1),''');']);
           while (Pr_temp ~= 1)&&(Pr_temp ~= 2)&&(Pr_temp ~= 3)         %检测输入的优先级是否在1/2/3范围内
                 disp(['请键入第',num2str(i),'个有效数据！优先级只能为1/2/3！']);
                 Pr_temp = input('Pr = ');
           end
           Pr(i) = Pr_temp;
       end
       Pr 
%{
       for i = 1:Tao
           disp(['请从τ1到τn依次键入各数据流的优先级Pr',num2str(i),'/共',num2str(Tao),'个 (1/2/3)']);
           Pr_temp = input('Pr = ');
           while (Pr_temp ~= 1)&&(Pr_temp ~= 2)&&(Pr_temp ~= 3)         %检测输入的优先级是否在1/2/3范围内
                 disp('请键入有效数据！优先级只能为1/2/3！');
                 Pr_temp = input('Pr = ');
           end
           %Pr(i) = input('Pr = ');
           Pr(i) = Pr_temp;
           disp(' ');
       end
%}
if ismember(3,Pr)       %如果有Class A的数据流
    %{
   以下代码被find函数代替
    A = 1;
    for i=1:length(Pr)
        if Pr(i) == 3       %Class A优先级最高是3
            j_a(A) = i;     %挑出Pr中的Class A
        end
            A = A + 1;
    end
    A = 1;
    for i=1:length(j_a)
        if j_a(i) ~= 0
            temp = j_a(i);
            J_A(A) = temp;      %挑出j_a中的非0元素，J_A为classA的所有数据流 
            A = A + 1;
        end
    end
    
    A = 1;
    for i = 1:length(j_a)
        if (j_a(i) ~= 0) && (j_a(i) ~= I)       %挑出J_A中的非I干扰流
            temp = j_a(i);
            j_Reala(A) = temp;
        end
        A = A + 1;
    end
    A = 1;
    for i=1:length(j_Reala)
        if j_Reala(i) ~= 0
            temp = j_Reala(i);
            J_RealA(A) = temp;      %挑出j_Reala中的非0元素，J_RealA为classA的除被研究对象的所有数据流  //[4]
            A = A + 1;              %这一步如果放到了if外面就不会起到除去0元素的作用
        end
    end
    %}
    J_A = find(Pr == 3);
    J_A_temp = find(Pr == 3);
    %if ismember(I,J_A)      %如果I是Class A
    J_A_temp(find(J_A_temp == I)) = [];   %删除I所在的位置的数据
        %location_I = find(J_A == I);
        %J_A(location_I) = [];
    J_RealA = J_A_temp;      %J_RealA必须存在
    %end
    clear location_I J_A_temp;
end

if ismember(2,Pr)       %如果有Class B的数据流
    %{
    以下代码被find函数代替
    B = 1;
    for i=1:length(Pr)
        if Pr(i) == 2       %Class B优先级是2
            j_b(B) = i;     %挑出Pr中的Class B
        end
        B = B + 1;
    end
    B = 1;
    for i=1:length(j_b)
        if j_b(i) ~= 0
            temp = j_b(i);
            J_B(B) = temp;      %挑出j_b中的非0元素，J_B为classB的所有数据流  //[1 2 3]
            B = B + 1;
        end
    end
    
    B = 1;
    for i = 1:length(j_b)
        if (j_b(i) ~= 0) && (j_b(i) ~= I)
            temp = j_b(i);
            j_Realb(B) = temp;  %去掉了与干扰流相同Class的被研究对象I
        end
        B = B + 1;
    end
    B = 1;
    for i=1:length(j_Realb)
        if j_Realb(i) ~= 0
            temp = j_Realb(i);
            J_RealB(B) = temp;      %挑出j_Realb中的非0元素，J_RealB为classB的除被研究对象的所有数据流  //[2 3]
            B = B + 1;              %这一步如果放到了if外面就不会起到除去0元素的作用
        end
    end
    %}
    J_B = find(Pr == 2);
    J_B_temp = find(Pr == 2);
    %if ismember(I,J_B)      %如果I是Class B
    J_B_temp(find(J_B_temp == I)) = [];
        %location_I = find(J_B == I);
        %J_B(location_I) = [];
    J_RealB = J_B_temp;
    %end
    clear location_I J_B_temp;
end
       
disp('以下数据直接键入。');
R = xlsread('data',1,'B2');    % Mbit/s
%R = input('\n传输速率:\nR=');
Alfa = R * [0.5 0.5 0.25 0.75];
global alfa_neg_A;
global alfa_pos_A;
global alfa_neg_B;
global alfa_pos_B;
       alfa_pos_A = Alfa(1);
       alfa_neg_A = -Alfa(2);
       alfa_pos_B = Alfa(3);
       alfa_neg_B = -Alfa(4);

      
%I = input('\n请键入被研究对象的角标\nI=');
global sl;
sl = xlsread('data',1,'D2');       
%sl = input('\n请键入switch latency\nsl = ');
global t;
t = xlsread('data',1,'E2');
%t=input('\n请键入I生成时间t\n');

%I与自己无法产生竞争
S_first_ii_max_i = 0;
S_last_ii_min_i = 0;
A_ii = 0;
%*************************************************************


%% ************【把S_In弄成三维矩阵并得到S_InSeq】*****************
%每一行代表一个输入端口，每一列没有实际含义（一行中的某一列代表输入端口中的一个帧或流），每一页代表一个交换机
%{
Num = input('\n交换机个数：\nNum = ');
Max_Input = input('\n最多几个输入端口：\nMax_Input = ');
Max_InputNum = input('\n某个输入端口最多的输入数量\nMax_InputNum = ');      %决定矩阵每一行元素数量
%}

%Num = input('\n交换机个数：\nNum = ');
Num = xlsread('data',1,'I2');

%Max_Input = input('\n最多几个输入端口：\nMax_Input = ');
Max_Input = xlsread('data',1,'J2');

%Max_InputNum = input('\n某个输入端口的，最多输入数量\nMax_InputNum = ');      %决定矩阵每一行元素数量
Max_InputNum = xlsread('data',1,'K2');
disp(['下列矩阵元素个数为',num2str(Max_InputNum),'个']);

global S_InSeq;               %设置输入序列
S_In = zeros(Max_Input,Max_InputNum,Num);  %行Max_Input代表最多几个输入端口，列Max_InputNum代表某个输入端口最多的输入数量，页Num代表交换机数量
%{
%判断是否需要删除旧的In_test.xls而生成新表
judge = exist('In_test.xls','file');
if judge == 2   %如果目前存在已有In_test.xls文件
    select = input('是否需要删除已有In_test.xls并重新生成？ Y/N\n','s');
    if select == 'Y'
        delete('In_test.xls');   %删除In.xls以方便后续重新生成最新In.xls表格
    end
end

%生成In_test.xls表格以键入数据
for i = 1:Num
    eval(['xlswrite(''In_test'',S_In(:,:,',num2str(i),'),',num2str(i),');']);    %Sheet i就代表Si的输入数据.把纯0矩阵写到In.xls中
    switch_loca = {'交换机位置'};
    eval(['xlswrite(''In_test'',switch_loca,',num2str(i),',''A',num2str(Max_Input + 2),''');']);    %把下面隔两行的位置处，设置输入交换机位置的格子
end

disp('请打开In_test.xls文件,并在Sheet X中输入SX的数据');
pause;
disp('所有交换机数据数据键入完成，666');
%}

%把In_test.xls表格内的内容读到S_In内。Sheet1为S1，Sheet2为S2
%旧版：
%{
 以下为旧版参数输入代码 
       for i = 1:Num
           disp('********************');
           for j = 1:Max_Input
               disp(['请键入S',num2str(i),'交换机的第',num2str(j),'个输入端口的数据，']);
               disp(['数据键入格式为行向量，保证“元素个数等于最多输入数量”，如：[a b c]。']);
               S_In(j,:,i) = input(' ');
               disp(' ');
           end
       end
%}
%新版：
for i = 1:Num
    for j = 1:Max_Input
            range = [num2str(j),':',num2str(j)];    %读取第j行
            S_In(j,:,i) = xlsread('In_test.xls',i,range);     %按一行一行读,得到真正S_In
    end
end
%S1=[1 0 0;2 3 6;0 0 0];
%S2=[1 2 3;4 0 0;5 0 0];
%S3=[2 0 0;3 0 0;6 0 0];

%实现：第1页是S3，第2页是S1，第3页是S2。得到S_InSeq
%旧版：
%{
旧版的整理顺序
for i = 1:Num
    disp(['这是第',num2str(i),'个交换机，请确认这是S几，并将这个数字代替x：S_In(:,:,x)']);
    x=input('');
    temp = S_In(:,:,x);
    S_InSeq(:,:,i) = temp;       %S_InSeq是具有拓扑顺序的交换机结构，其第一页就代表第一个交换机
    disp(' ');
end
%}
%新版：
for i = 1:Num
    range = [num2str(1),':',num2str(Max_Input)];
    temp(:,:) = xlsread('In_test',i,range);     %读取了S1、S2、S3的输入
    seq = xlsread('In_test',i,['A',num2str(Max_Input + 3)]);     %seq代表这个交换机的位置（第几个）
    S_SwitchSeq(i) = seq;   %弄出S几是第几个交换机。
    S_InSeq(:,:,seq) = temp;
end
%S_InSeq(:,:,1) = S_In(:,:,3);
%S_InSeq(:,:,2) = S_In(:,:,1);
%S_InSeq(:,:,3) = S_In(:,:,2);

for i = 1:Num
    n = 1;
    while ~isequal(S_In(:,:,n),S_InSeq(:,:,i))
        n = n + 1;
    end
    S_ReSwitchSeq(i) = n;   %弄出第几个交换机是S几。
end

clear judge select switch_loca title1 title2 range temp seq i j n;


%% **************【整理S_Out与S_OutSeq】***********************
%{
Max_Output = input('\n最多几个输出端口：\nMax_Output = ');
Max_OutputNum = input('\n某个输出端口最多的输出数量\nMax_OutputNum = ');     %决定矩阵每一行元素数量     
S_Out = zeros(Max_Output,Max_OutputNum,Num);  %行代表最多几个输出端口，列代表某个输出端口最多的输出数量，页代表交换机数量
       for i = 1:Num        %为了得到每一个交换机的输出序列
           for j = 1:Max_Output
               disp(' ');
               disp(['请键入S',num2str(i),'交换机的输出数据，']);
               disp('数据键入格式为行向量，保证“元素个数等于最多输出数量”，如：[a b c]。');
               S_Out(j,:,i) = input(' ');
               disp(' ');
           end
       end
%}
%Max_Output = input('\n最多几个输出端口：\nMax_Output = ');
Max_Output = xlsread('data',1,'L2');

%Max_OutputNum = input('\n某个输出端口的，最多输出数量\nMax_OutputNum = ');     %决定矩阵每一行元素数量
Max_OutputNum = xlsread('data',1,'M2');
disp(['下列矩阵元素个数为',num2str(Max_OutputNum),'个']);

global S_OutSeq;               %设置输出序列
S_Out(:,:,1) = [1 2 3 0 0;6 0 0 0 0];
S_Out(:,:,2) = [1 2 3 4 5;0 0 0 0 0];
S_Out(:,:,3) = [2 3 6 0 0;0 0 0 0 0];
%{
%实现：第1页是S3，第2页是S1，第3页是S2。
for i = 1:Num
    disp(['这是第',num2str(i),'个交换机，请确认这是S几，并将这个数字代替x：S_Out(:,:,x)']);
    x=input('');
    temp = S_Out(:,:,x);
    S_OutSeq(:,:,i) = temp;       %S_InSeq是具有拓扑顺序的交换机结构，其第一页就代表第一个交换机
    S_OutSeqNo(i) = x;
    disp(' ');
end
%}
S_OutSeq(:,:,1) = S_Out(:,:,3);
S_OutSeq(:,:,2) = S_Out(:,:,1);
S_OutSeq(:,:,3) = S_Out(:,:,2);
S_OutSeqNo(1) = 3;
S_OutSeqNo(2) = 1;
S_OutSeqNo(3) = 2;
%不清除S_Out是因为要留下直接用S几来确定输入输出参数


%% ***************【寻找Pi，即I的路径】************************
%路径矩阵的第一行代表交换机序号S几，每一列代表一个交换机上的输出端口
p = 0;
for i = 1:Num
    for j = 1:Max_Output
        m = 0;
        for k = 1:Max_OutputNum
            if S_OutSeq(j,k,i) == I     %有被研究对象存在，意味着这个输出端口属于Pi
                m=m+1;
                p=p+1;
            end
            if m == 1
                Pi(:,p) = [S_OutSeqNo(i),j];        %Pi为路径矩阵（第一行是交换机序号）
            end
        end
    end
end
clear i j m p;


%% **********************【找出source_i/j】************************************
%这是找到每一个数据流的出生地位置（以交换机为起点，比如2数据流出生在S3，第一个交换机，记录的就是1）
source = zeros(Tao,Num);    %一共有Tao个数据流，最多Num个有输入的交换机。这是为了防止：
                            %从上一级交换机只传递了一个数据流到下一交换机
count_source = 1;
for i = 1:Num
    for j = 1:Max_Input
        n = 0;
        for k = 1:Max_InputNum
            if S_In(j,k,i) ~= 0
                n = n + 1;
            end
        end
        if n == 1    %这一个输入端口中只有一个输入。
            loca_source = find(S_In(j,:,i) ~= 0);    %找到每一个输入端口中的单独输入的位置（避免人为输入的随意性）
            source_flow = S_In(j,loca_source,i);     %输入端口中的单独输入的角标
            count_source = 1;
            while ~isequal(S_In(:,:,i),S_InSeq(:,:,count_source))     
                count_source = count_source + 1;      
            end
            source(source_flow,i) = count_source;     %source矩阵意义：
                                                      %数值代表第n个交换机，列数代表交换机序号，
                                                      %行数代表数据流角标
        end
    end
end
for i = 1:Tao
    temp = find(source(i,:) ~= 0);
    Flow_Source(i) = source(i,temp);    %Flow_Source代表了每一个数据流的起始位置
end
clear source source_flow loca_source count_source i temp;


%% **********************【找出first_ij】*************************************
%把Pi上的所有输出扫描一遍，
%把具有i和j的交换机在S_InSeq中排个序，把第一个拿出来就是first_ij，最后一个就是last_ij
disp(' ');
disp('****************************');
disp('【请记录】：first_ij与last_ij：');
%找到所有Class B的数据流
if ismember(2,Pr)
    judge_first = 0;                %判断是否为first_ij的判断参数
    for J = 1:length(J_B)
        for i = 1:length(Pi)        %从Pi的第一个交换机到Pi的最后一个交换机
            %需要弄出i=1时是哪个交换机。
            Switch_No = Pi(1,i);    %Pi上的第i个交换机是S【Switch_No】，也就是S_Out中的第【Switch_No】页
            for j = 1:Max_Output
                m = 0;
                %{
                不用这个功能是因为无法排除掉first_11【I = J_B(J)】这个因素。
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
                if m == 2                % m==2说明这个输出中既有I，也有J_B(J),即循环轮到的干扰流
                    %cmpt(n) = i;         %(cmpt是compete的意思，代表第几个交换机)
                    %cmpt_out(n) = j;     %(cmpt_out代表某个交换机上，竞争存在的第几个输出端口)
                    %n = n + 1;
                    m = 0;
                    judge_first = judge_first + 1;      %judge_first自增1
                end
            end
            if judge_first == 1;                        %在第一次增加的时候，也就是第一次输出中既有I也有J_B的时候，是first_ij
                first_ij = Switch_No;                   %代表交换机S的序号
                for l = 1:Max_Output
                    if ismember(I,S_Out(l,:,first_ij))
                        first_ij_port = l;
                    end
                end
                disp(['first_',num2str(I),num2str(J_B(J)),' = S',num2str(first_ij),num2str(first_ij_port),'交换机是S',num2str(first_ij)]);
                %利用结构体Cmpt.first来存储first_ij
                %查询first_ij
                %不考虑I = J_B(J)的情况
                Cmpt.first(I,J_B(J),1) = first_ij;  %行代表I，列代表j，第一页代表交换机序号，第二页代表端口号
                Cmpt.first(I,J_B(J),2) = first_ij_port;
                break
            end
        end
        judge_first = 0;        %把Pi走完一趟之后judge_first清零，换下一个优先级的干扰
    end
    clear i j k m l first_ij_port;
end

%找所有的Class A流
if ismember(3,Pr)
    judge_first = 0;                %判断是否为first_ij的判断参数
    for J = 1:length(J_A)
        [hang lie] = size(Pi);
        for i = 1:lie        %从Pi的第一列到Pi的最后一列
            %需要弄出i=1时是哪个交换机。
            Switch_No = Pi(1,i);    %Pi上的第i个交换机是S“Switch_No”，也就是S_Out中的第【Switch_No】页
            for j = 1:Max_Output
                m = 0;
                if ismember(J_A(J),S_Out(j,:,Switch_No))
                    for k = 1:Max_OutputNum
                        if (S_Out(j,k,Switch_No) == I) || (S_Out(j,k,Switch_No) == J_A(J))
                            m=m+1;
                        end
                    end
                end
                if m == 2                % m==2说明这个输出中既有I，也有J_A(J),即循环轮到的干扰流
                    %cmpt(n) = i;         %(cmpt是compete的意思，代表第几个交换机)
                    %cmpt_out(n) = j;     %(cmpt_out代表某个交换机上，竞争存在的第几个输出端口)
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
                disp(['first_',num2str(I),num2str(J_A(J)),' = S',num2str(first_ij),num2str(first_ij_port),'交换机是S',num2str(first_ij)]);
                %利用结构体Cmpt.first来存储first_ij
                Cmpt.first(I,J_A(J),1) = first_ij;  %行代表I，列代表j，第一页代表交换机序号，第二页代表端口号
                Cmpt.first(I,J_A(J),2) = first_ij_port;
                break
            end
        end
        judge_first = 0;
    end
    clear i j k m J judge_first first_ij Switch_No l first_ij_port hang lie;
end


%% *******************【找出last_ij】**************************
%{
本来想通过在某个口同时存在I&J，由此进入监视状态，
若二者继续同时存在，则判断为真，通过；
若条件消失，则上一个交换机为last，
然后再测一次是在上一个交换机的哪个端口。
但是监视状态的进入和退出很难实现！

新思路：直接逆序寻找first就好了。
%}
reverse_Pi = fliplr(Pi);
%找到所有Class B的数据流
if ismember(2,Pr)
    judge_last = 0;                         %判断是否为last_ij的判断参数
    for J = 1:length(J_B)
        for i = 1:length(reverse_Pi)
            Switch_No = reverse_Pi(1,i);    %reverse_Pi上的第i个交换机是S“Switch_No”，也就是S_Out中的第【Switch_No】页
            for j = 1:Max_Output    %端口数量
                m = 0;
                if ismember(J_B(J),S_Out(j,:,Switch_No))
                    for k = 1:Max_OutputNum
                        if (S_Out(j,k,Switch_No) == I) || (S_Out(j,k,Switch_No) == J_B(J))
                            m=m+1;
                        end
                    end
                end
                if m == 2      % m==2说明这个输出中既有I，也有J_B(J),即循环轮到的干扰流
                    %{
                    %【开始进入last_ij的监视状态，若继续二者同时存在则判断为真——通过；】
                    %【若不同时存在——则上一个交换机为last，再测一次在上一个交换机的哪个端口】
                    whatchdog = 1;      %启动监视
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
                disp(['last_',num2str(I),num2str(J_B(J)),' = S',num2str(last_ij),num2str(last_ij_port),'交换机是S',num2str(last_ij)]);
                %利用结构体Cmpt.last来存储last_ij
                %查询last_ij
                %不考虑I = J_B(J)的情况
                Cmpt.last(I,J_B(J),1) = last_ij;
                Cmpt.last(I,J_B(J),2) = last_ij_port;
                break
            end
        end
        judge_last = 0;
    end
    clear J i j m k l judge_last last_ij_port;
end

%找到所有Class A的数据流
if ismember(3,Pr)
    judge_last = 0;                         %判断是否为last_ij的判断参数
    for J = 1:length(J_A)
        [hang lie] = size(reverse_Pi);
        for i = 1:lie
            Switch_No = reverse_Pi(1,i);    %Pi上的第i个交换机是S“Switch_No”，也就是S_Out中的第【Switch_No】页
            for j = 1:Max_Output    %端口数量
                m = 0;
                if ismember(J_A(J),S_Out(j,:,Switch_No))
                    for k = 1:Max_OutputNum
                        if (S_Out(j,k,Switch_No) == I) || (S_Out(j,k,Switch_No) == J_A(J))
                            m=m+1;
                        end
                    end
                end
                if m == 2     % m==2说明这个输出中既有I，也有J_A(J),即循环轮到的干扰流
                    %{
                    %【开始进入last_ij的监视状态，若继续二者同时存在则判断为真——通过；】
                    %【若不同时存在——则上一个交换机为last，再测一次在上一个交换机的哪个端口】
                    whatchdog = 1;      %启动监视
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
                disp(['last_',num2str(I),num2str(J_A(J)),' = S',num2str(last_ij),num2str(last_ij_port),'交换机是S',num2str(last_ij)]);
                %利用结构体Cmpt.last来存储last_ij
                %查询last_ij
                %不考虑I = J_A(J)的情况
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


%% ************************************************************
%{
%Cmpt(1,:) = cmpt;
%Cmpt(2,:) = cmpt_out;
%Cmpt
%手动对其进行排序。
%CMPT = input('\n请参考拓扑结构图，根据存在I、J竞争的交换机的顺序，对比Cmpt的第一行元素序号，以矩阵的形式，以一列为单位，按先后顺序键入新矩阵\n'); %Cmpt第一行的数字是几就代表S几
%CMPT = [1 2;1 1];
%CMPT1 = size(CMPT);  %取CMPT的规格，行数固定为2，列数代表存在竞争的输出端口数

disp(' ');
disp('****************************');
disp(['first_ij = S',num2str(CMPT(1,1)),num2str(CMPT(2,1))]);       %CMPT(1,1)代表交换机S几，CMPT(2,1)代表输出端口几
disp(['last_ij = S',num2str(CMPT(1,CMPT1(2))),num2str(CMPT(2,CMPT1(2)))]);
disp('请记录下first_ij与last_ij');
disp('****************************');
disp(' ');
%}


%% ************************不知道有没有用***********************
%{
不知道有没有用
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


%% ********************* 计算【S_first_ij_min_j】和【S_first_ij_min_i】的值【还缺针对A的计算】 **************************
%{
%先显示first_ij是多少，如S11
%disp(['first_ij = S',num2str(CMPT(1,1)),num2str(CMPT(2,1))]);       %CMPT(1,1)代表交换机S几，CMPT(2,1)代表输出端口几
%再搞出first_ij的交换机是第No_first_ij个，如第“2”个，要找到CMPT(1,1)的位置
%n = 1;
%while ~isequal(S_In(:,:,CMPT(1,1)),S_InSeq(:,:,n))      %看S【CMPT(1,1)】与拓扑结构中第几个交换机相等
%        n = n + 1;      %当找到是first_ij是第n个交换机的时候，跳出循环。
%end
%No_first_ij = n;
%}

%S.max = zeros(Max_Output,Tao,Num);   %保存所有的S参数，页——I/J,数量上 = τ，行——交换机名称，
%S.min = zeros(Max_Output,Tao,Num);   %列——输出端口数量。
%对于ClassB，搞出J从第No_J个交换机进去的，如第“1”个【先找J，再找非J元素】
%***********************************************************
disp('这是针对J属于Class B的计算：');
count_j = 1;
for J=1:length(J_RealB)
    for i=1:Num                             %页
        for j=1:Max_Input                   %行
            n = 0;
            N = 0;
            for k=1:Max_InputNum            %列
                if S_InSeq(j,k,i) == J_RealB(J)      %一个一个找J_B
                    n=n+1;
                end
                if S_InSeq(j,k,i) == I      %一个一个找I
                    N=N+1;
                end
            end
            %********************
            if n == 1                       %代表这一行（这个输入端口）有J存在
                m = 0;
                for l = 1:Max_InputNum
                    if (S_InSeq(j,l,i) ~= J_RealB(J)) && (S_InSeq(j,l,i) ~= 0)     %在这一行（这个输入端口）有非J数据流存在
                        m=m+1;
                    end
                end
                if m == 0                   %不存在非J元素（只有J存在）
                    No_J_B(count_j) = i;      %No_J这个数组表示S_InSeq中的第几个交换机
                    count_j = count_j + 1;
                end
            end
            %********************
            if N == 1                       %代表这一行（这个输入端口）有I存在
                M = 0;
                for L=1:Max_InputNum
                    if (S_InSeq(j,L,i) ~= I) && (S_InSeq(j,L,i) ~= 0)     %在这一行（这个输入端口）有非I数据流存在
                        M=M+1;
                    end

                end
                if M == 0                %不存在非I元素（只有I存在）
                    No_I = i;
                end
            end
            %********************
        end
    end
end

%计算S_first_ij_min_j
%S_first_ij_min_j = C（I）*（No_first_ij - No_J_B + 1）
disp('**************【计算S_first_ij_min_j】*******************');
disp('********************************************************');
for count_j = 1:length(J_RealB)
    %Name_first_ij = input(['请按顺序键入之前记录的first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'的交换机序号和端口，以行向量的形式键入：\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealB(count_j),1)),S_InSeq(:,:,n))      %看S【CMPT(1,1)】与拓扑结构中第几个交换机相等
        n = n + 1;      %当找到是first_ij是第n个交换机的时候，跳出循环。
    end
    %存储first_ij的位置、序号（名字）、端口
    %见S.min的参数
    No_first_B(I,J_RealB(count_j)) = n;   
    Name_first_B(I,J_RealB(count_j)) = Cmpt.first(I,J_RealB(count_j),1);
    Port_first_B(I,J_RealB(count_j)) = Cmpt.first(I,J_RealB(count_j),2);
    %No_first_ij = input(['first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'所在的S',...
    %    num2str(Cmpt.first(I,J_RealB(count_j),1)),'是第几个交换机：\nNo_first_ij = ']);
    %disp('【请记录：】');
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

%计算S_first_ij_min_i
%S_first_ij_min_i = C(I) * (No_first_ij - No_I + 1)
disp('**************【计算S_first_ij_min_i】*******************');
disp('********************************************************');
disp(['这是I的计算，I = ',num2str(I)]);
if Pr(I) == 3
    disp('I是Class A的数据流');
else if Pr(I) == 2
        disp('I是Class B的数据流');
    else
        disp('I是Class C的数据流');
    end
end
disp(' ');
for count_j = 1:length(J_RealB)
    %Name_first_ij = input(['请按顺序键入之前记录的first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'的交换机序号和端口，以行向量的形式键入：\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealB(count_j),1)),S_InSeq(:,:,n))      
        n = n + 1;      %当找到是first_ij是第n个交换机的时候，跳出循环。
    end
    No_first_B(I,J_RealB(count_j)) = n;
    %No_first_ij = input(['first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'所在的S',...
    %    num2str(Cmpt.first(I,J_RealB(count_j),1)),'是第几个交换机：\nNo_first_ij = ']);
    %disp('【请记录：】');
    %disp('--------------------------');
    %eval(['S_first_',num2str(I),num2str(J_RealB(count_j)),'_min_',num2str(I),' = ',...
    %    num2str(C(J_RealB(count_j)) * (No_first(I,J_RealB(count_j)) - No_I + 1))]);
    disp('--------------------------');
    disp(' ');
    %S.min(Name_first_ij(1),Name_first_ij(2),I) =...
    %    C(J_RealB(count_j)) * (No_first_ij - No_I + 1);
    %S.min( Cmpt.first(I,J_RealB(count_j),1) , Cmpt.first(I,J_RealB(count_j),2) , I) =...
    %    C(J_RealB(count_j)) * (No_first_ij - No_I + 1);     %S.min是以S_mn为单位的存储空间，而不是以first/last为单位的存储空间
    S.min( Cmpt.first(I,J_RealB(count_j),1) , Cmpt.first(I,J_RealB(count_j),2) , I) =...
        C(J_RealB(count_j)) * (No_first_B(I,J_RealB(count_j)) - No_I + 1);     %S.min是以S_mn为单位的存储空间，而不是以first/last为单位的存储空间
end
disp('*********************************');
disp(' ');
disp(' ');

%***********************************************************
disp('这是针对J属于Class A的计算：');
count_j = 1;
for J=1:length(J_RealA)
    for i=1:Num                             %页
        for j=1:Max_Input                   %行
            n = 0;
            N = 0;
            for k=1:Max_InputNum            %列
                if S_InSeq(j,k,i) == J_RealA(J)      %一个一个找J_A
                    n=n+1;
                end
                if S_InSeq(j,k,i) == I      %一个一个找I
                    N=N+1;
                end
            end
            %********************
            if n == 1                       %代表这一行（这个输入端口）有J存在
                m = 0;
                for l = 1:Max_InputNum
                    if (S_InSeq(j,l,i) ~= J_RealA(J)) && (S_InSeq(j,l,i) ~= 0)     %在这一行（这个输入端口）有非J数据流存在
                        m=m+1;
                    end
                end
                if m == 0                   %不存在非J元素（只有J存在）
                    No_J_A(count_j) = i;      %No_J这个数组表示S_InSeq中的第几个交换机
                    count_j = count_j + 1;
                end
            end
            %********************
            if N == 1                       %代表这一行（这个输入端口）有I存在
                M = 0;
                for L=1:Max_InputNum
                    if (S_InSeq(j,L,i) ~= I) && (S_InSeq(j,L,i) ~= 0)     %在这一行（这个输入端口）有非I数据流存在
                        M=M+1;
                    end

                end
                if M == 0                %不存在非I元素（只有I存在）
                    No_I = i;
                end
            end
            %********************
        end
    end
end

%计算S_first_ij_min_j
%S_first_ij_min_j = C（I）*（No_first_ij - No_J_A + 1）
disp('**************【计算S_first_ij_min_j】*******************');
disp('********************************************************');
for count_j = 1:length(J_RealA)
    %Name_first_ij = input(['请按顺序键入之前记录的first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'的交换机序号和端口，以行向量的形式键入：\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealA(count_j),1)),S_InSeq(:,:,n))      %对比S_In与S_InSeq
        n = n + 1;      %当找到是first_ij是第n个交换机的时候，跳出循环。
    end
    %存储first_ij的位置、序号（名字）、端口
    %见S.min的参数
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

%计算S_first_ij_min_i
%S_first_ij_min_i = C(I) * (No_first_ij - No_I + 1)
disp('**************【计算S_first_ij_min_i】*******************');
disp('********************************************************');
disp(['这是I的计算，I = ',num2str(I)]);
if Pr(I) == 3
    disp('I是Class A的数据流');
else if Pr(I) == 2
        disp('I是Class B的数据流');
    else
        disp('I是Class C的数据流');
    end
end
disp(' ');
for count_j = 1:length(J_RealA)
    %Name_first_ij = input(['请按顺序键入之前记录的first_',num2str(I),...
    %    num2str(J_RealA(count_j)),'的交换机序号和端口，以行向量的形式键入：\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealA(count_j),1)),S_InSeq(:,:,n))     
        n = n + 1;      %当找到是first_ij是第n个交换机的时候，跳出循环。
    end
    No_first_A(I,J_RealA(count_j)) = n;
    disp('--------------------------');
    disp(' ');
    S.min( Cmpt.first(I,J_RealA(count_j),1) , Cmpt.first(I,J_RealA(count_j),2) , I) =...
        C(J_RealA(count_j)) * (No_first_A(I,J_RealA(count_j)) - No_I + 1);     %S.min是以S_mn为单位的存储空间，而不是以first/last为单位的存储空间
end
disp('*********************************');
disp(' ');
disp(' ');
clear i j k l L m M n N count_j J;


%% ***********************【计算S_first_ij_max_j和S_first_ij_max_i】*******************
 %{

%初步估计：应该沿着Pi计算
CreditA = 0;
CreditB = 0;
[hang lie] = size(Pi);  %hang没有什么含义，lie代表一个交换机上的输出端口（先决定lie，再决定hang）
%for i = 1:lie
    i = 1;
    S_research_temp = S_Out(Pi(2,i),:,Pi(1,i));  %对于S几这个交换机输出端口的输出数据，这是个一维数组  [2 0 3 6 0]
    S_research_local = find(S_research_temp ~= 0);  %对于某个输出端口来说S_research_local = [1 3 4]
    for localS = 1:length(S_research_local)         %
        %S_research = [2 3 6]  无0元素的输出数据
        
        S_research(localS) = S_Out(Pi(2,i),S_research_local(localS),Pi(1,i));   %localS = 2时，S_research_local = 3
    end

        %clear S_research_temp S_research_local localS
        %*****************************
        %开始分Class
        ClassA = [];
        ClassB = [];
        ClassC = [];
        for k = 1:length(S_research)
            Pr_research(k) = Pr(S_research(k));
        end
        %for k = 1:length(S_research)
        %    Pr_research(k) = Pr(S_research(k)); %获得该输出端口所有数据的优先级
            ClassA_temp = find(Pr_research == 3);   %获得该输出端口中优先级=3的角标位置
            for localA = 1:length(ClassA_temp)
            ClassA(localA) = S_Out(Pi(2,i),ClassA_temp(localA),Pi(1,i));    %获得ClassA的数组
            end
        %end
        %for k = 1:length(S_research)
        %    Pr_research(k) = Pr(S_research(k)); %获得该输出端口所有数据的优先级
            ClassB_temp = find(Pr_research == 2);   %获得该输出端口中优先级=2的角标位置
            for localB = 1:length(ClassB_temp)
                ClassB(localB) = S_Out(Pi(2,i),ClassB_temp(localB),Pi(1,i));
            end
        %end
        %for k = 1:length(S_research)
        %    Pr_research(k) = Pr(S_research(k)); %获得该输出端口所有数据的优先级
            ClassC_temp = find(Pr_research == 1);   %获得该输出端口中优先级=1的角标位置
            for localC = 1:length(ClassC_temp)
                ClassC(localC) = S_Out(Pi(2,i),ClassC_temp(localC),Pi(1,i));
            end
        %end
        %clear ClassA_temp ClassB_temp ClassC_temp localA localB localC k
        %*****************************
        %把S_research中的第j个元素放到其所在Class数组中的最后一个
    for j = 1:length(S_research)
        delaytime = 0;
        if Pr(S_research(j)) == 3
            a = find(ClassA == S_research(j));  %S_research中的第j个元素在ClassA中的位置
            temp = ClassA(a);
            ClassA(a) = []; %删除这个元素
            ClassA = [ClassA temp]; %将这个第j个元素置到最后
        else if Pr(S_research(j)) == 2
                b = find(ClassB == S_research(j));  %S_research中的第j个元素在ClassA中的位置
                temp = ClassB(b);
                ClassB(b) = []; %删除这个元素
                ClassB = [ClassB temp]; %将这个第j个元素置到最后
             else if Pr(S_research(j)) == 1
                c = find(ClassC == S_research(j));  %S_research中的第j个元素在ClassA中的位置
                temp = ClassC(c);
                ClassC(c) = []; %删除这个元素
                ClassC = [ClassC temp]; %将这个第j个元素置到最后
                end
            end
        end
        %clear a b c temp
        %*****************************
        %取刚才三个Class数组各自的C和Pr
        for tiqu = 1:length(ClassA) %tiqu意味：提取
            Pr_A(tiqu) = Pr(ClassA(tiqu));
            C_A(tiqu) = C(ClassA(tiqu));
        end
        for tiqu = 1:length(ClassB)
            Pr_B(tiqu) = Pr(ClassB(tiqu));
            C_B(tiqu) = C(ClassB(tiqu));
        end
        for tiqu = 1:length(ClassC)
            Pr_C(tiqu) = Pr(ClassC(tiqu));
            C_C(tiqu) = C(ClassC(tiqu));
        end
        %*****************************
        %正式开始计算
        %目前只计算只有B的情况
        if ~isempty(ClassA)     %如果ClassA数组不是空集
            disp('暂时不考虑这种情况');   
        else if ~isempty(ClassB)        %如果ClassA数组是空集 且 ClassB数组不是空集
                cc = 1;     %在ClassC中计数
                for trans = 1:(length(ClassB)-1)
                    if CreditB>=0
                        %{
                        if ~isempty(ClassC) && trans == length(ClassB)  %从被研究对象数组中的倒数插空
                            delaytime = delaytime + C(ClassC(cc));
                            cc = cc + 1;
                            CreditB = CreditB + C(ClassC(cc)) * alfa_pos_B;
                        end
                        %}
                        delaytime = delaytime + C(ClassB(trans));   %传输ClassB的帧
                        CreditB = CreditB + C(ClassB(trans)) * alfa_neg_B; %传输完成后的CreditB
                    end
                    if CreditB<0
                        Replenish_B = (0 - CreditB) / (alfa_pos_B); %从负值恢复到0的用时
                        if ~isempty(ClassC) %&& (Replenish_B >= C(ClassC(cc))) 如果C非空

                            delaytime = delaytime + Replenish_B;    %CreditB在恢复同时，并不传输C
                            CreditB = CreditB + Replenish_B * alfa_pos_B;
                            cc = cc + 1;
                        else if ~isempty(ClassC) && (Replenish_B < C(ClassC(cc)))
                                delaytime = delaytime + C(ClassC(cc));  % non-preemptive
                                CreditB = CreditB + C(ClassC(cc)) * alfa_pos_B;
                                cc = cc + 1;
                            else if isempty(ClassC)
                                    delaytime = delaytime + Replenish_B;
                                    CreditB = CreditB + Replenish_B * alfa_pos_B;
                                end
                            end
                        end
                    end
                    
                end
            end
        end
        S.max(Pi(1,i),Pi(2,i),j) = delaytime + C(ClassB(j));
        
        %现存问题：如何让1跟2、3不一样？！！！！！！！！！！！
    end
%end
a = 0;
if a>-1
disp(' ');
end
 %}


S.max(1,1,1) = 40;
S.max(1,1,2) = 360;
S.max(1,1,3) = 360;
S.max(2,1,4) = 40;


%{
%S.max(3,1,3) = 40;%不一定对（对I=1不存在这一行）
%{

%% ***********************【计算S_first_ij_max_j和S_first_ij_max_i】*******************

%S_S11_min_2
S_first_ij_min_j = C(J_b) * (No_first_ij - No_J + 1);
%S_S11_min_1
S_first_ij_min_i = C(I) * (No_first_ij - No_I + 1);

%***********************************************
%计算max有待解决（根据TA）
%S_S11_max_2
S_first_ij_max_j = 360;

%S_S11_max_1
S_first_ij_max_i = 40;
%***********************************************

%M_S11_1 = 40
%I在firsti，j的最早开始时间
M_first_ij_i = C(I) * (No_first_ij - No_I + 1);


%M_ECU1_1 = 0
M_first_i_i = 0;


%M_S21_1 = 80
%显示一下last_ij是多少
disp(['last_ij = S',num2str(CMPT(1,CMPT1(2))),num2str(CMPT(2,CMPT1(2)))]);
%再搞出last_ij的交换机是第No_last_ij个，如第“3”个，要找到“CMPT(1,CMPT1(2))”的位置
n = 1;
while ~isequal(S_In(:,:,CMPT(1,CMPT1(2))),S_InSeq(:,:,n))      %比较两个矩阵是否相等
    n = n + 1;         %当找到是last_ij是第n个交换机的时候，跳出循环。
end
No_last_ij = n;        %last_ij是第“No_last_ij”个交换机
%计算M_last_i_i
M_last_i_i = C(I) * (No_last_ij - No_I+ 1);

%% ① D_hpX_i_t
A_i_j = S_first_ij_min_j - M_first_ij_i;
if Pr(I) == 1
    X = 'Class A';
else if Pr(I) == 2
        X = 'Class B';
    else
        X = 'Class C';
    end
end
disp(['因为被研究对象角标i=',num2str(I),'，所以τ',num2str(I),'属于',X]);
D_classX_IJt = positive(1 + floor((t + S_first_ij_max_i - S_first_ij_min_j + A_i_j) / T(J_b))) * C(J_b);
%}

%}


%% ***********************【M_h_i】*************************
%【错了！计算M应该按照first_ij来算】
%{
%计算M应该是要按着Pi来算  <----  是错的想法
for count_m = 1:(length(Pi))     %到后面跟ECU的数据进行数组拼接
    %m_h_i是不计算ECU的最早开始时间
    m_h_i(I,count_m) = (S_SwitchSeq(Pi(1,count_m)) - Flow_Source(I) + 1) * C(I);     %第一页是ECU，第二页开始是第一个交换机，以此类推
end
%M_h_i是计算ECU的最早开始时间
%以下这几步是为了应对I等于"某大于1的数"从而导致的m_h_i矩阵维度变化问题
temp = size(m_h_i);
temp1 = zeros(temp(1),1);
temp2 = [temp1,m_h_i]; 
M_h_i = temp2(temp(1),:);
clear count_m m_h_i temp temp1 temp2;
%}

%计算与B竞争的M
for count_m = 1:length(J_RealB)
    %m_h_i是不计算ECU的最早开始时间
    %存储数据用M_first_ij_i的形式，调用数据用M_S_mn_i的形式
    M_h_i_B(Cmpt.first(I,J_RealB(count_m),1) , Cmpt.first(I,J_RealB(count_m),2) , I) =...
        (S_SwitchSeq(Cmpt.first(I,J_RealB(count_m),1)) - Flow_Source(I) + 1) * C(I);
end
%计算与A竞争的M
for count_m = 1:length(J_RealA)
    %m_h_i是不计算ECU的最早开始时间
    %存储数据用M_first_ij_i的形式，调用数据用M_S_mn_i的形式
    M_h_i_A(Cmpt.first(I,J_RealA(count_m),1) , Cmpt.first(I,J_RealA(count_m),2) , I) =...
        (S_SwitchSeq(Cmpt.first(I,J_RealA(count_m),1)) - Flow_Source(I) + 1) * C(I);
end
M_ECU_i = 0;

clear count_m 


%% **************************【A_ij】**************************
%计算Class B的A_ij
for i = 1:length(J_RealB)
    % A_ij = S_first_ij_max_j - M_first_ij_i
    % Name_first(I,J_RealB(i))代表了交换机序号（名称），也就是S_mn中的m
    % Port_first(I,J_RealB(i))代表了输出端口，也就是S_mn中的n
    A(I,J_RealB(i)) = S.max(Name_first_B(I,J_RealB(i)) , Port_first_B(I,J_RealB(i)) , J_RealB(i)) -...
        M_h_i_B(Name_first_B(I,J_RealB(i)) , Port_first_B(I,J_RealB(i)) , I);
end
%计算Class A的A_ij
for i = 1:length(J_RealA)
    A(I,J_RealA(i)) = S.max(Name_first_A(I,J_RealA(i)) , Port_first_A(I,J_RealA(i)) , J_RealA(i)) - ...
        M_h_i_A(Name_first_A(I,J_RealA(i)) , Port_first_A(I,J_RealA(i)) , I);
end


%% ***********************【D_ClassX】************************
%直接计算I所在的那一组
%将自己这一组排除
if ismember(I,J_B)
    for i = 1:length(J_RealB)
        %如果I是Class B，则将I从其中先排除。
        D_ClassB(I,J_RealB(i)) = positive(1 + floor((t + S.max(Name_first_B(I,J_RealB(i)) , Port_first_B(I,J_RealB(i)) , I) -...
            S.min(Name_first_B(I,J_RealB(i)) , Port_first_B(I,J_RealB(i)) , J_RealB(i)) + A(I,J_RealB(i)))...
            / T(J_RealB(i)))) * C(J_RealB(i));
    end
    %单独计算I自己
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
        %单独计算I自己
        D_ClassA_I = positive(1 + floor((t + S_first_ii_max_i - S_last_ii_min_i + A_ii)...
            / T(I))) * C(I);
        D_ClassA(I,I) = D_ClassA_I;
        Sum_D_ClassA = sum(sum(D_ClassA));
        clear D_ClassA_I i;
    end
end


%% *****************************【Rep_i_X】*******************************
if ismember(I,J_B)
    Rep_B(I) = Sum_D_ClassB * (alfa_neg_B / alfa_pos_B);
else if ismember(I,J_A)
        Rep_A(I) = Sum_D_ClassA * (alfa_neg_A / alfa_pos_A);
    end
end


%% ************************【含有迭代变量的参数】********************************
%以下对于hp而言，只有X == B的时候才成立。所以其实是在计算D_hp_B
if ismember(I,J_B)
    %要先弄出所有的last_ij（j属于hpX，即Class A。因为此处不考虑I是Class C）
    Temp = 0;   %此值最后为D_hp_B
    for i = 1:length(J_A)
        disp(['last_',num2str(I),num2str(J_A(i)) , '是S' , num2str(Cmpt.last(I,J_A(i),1)) ,...
            num2str(Cmpt.last(I,J_A(i),2))]);
        eval(['syms W_S',num2str(Cmpt.last(I,J_A(i),1)),num2str(Cmpt.last(I,J_A(i),2)),...
            '_',num2str(I)]);
        %不知道要怎么处理
        %eval(['temp = positive(1 + floor( (W_S',num2str(Cmpt.last(I,J_A(i),1)),num2str(Cmpt.last(I,J_A(i),2)),...
        %    '_',num2str(I),') / T(J_A(i)) )) * C(J_A(i))']);
    end
end


%% ***********************【子函数】********************************
%这部分是以上代码所使用过的子函数
%*********************************************************
% I - TRANST


%*********************************************************





%% ③ D_lpX_i


%% ④ max(Cj)


%% ⑤ (num(Pi)-1)*sl


%% ⑥ delta_h_i_t


%% ⑦ max(alfa_neg_X / alfa_pos_X)

%}