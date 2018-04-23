%现在要合并S_In与S_Out

%% S_In
%每一行代表一个输入端口，每一列没有实际含义（一行中的某一列代表输入端口中的一个帧或流），每一页代表一个交换机

%Num = input('\n交换机个数：\nNum = ');
Num = xlsread('data',1,'I2');

%Max_Input = input('\n最多几个输入端口：\nMax_Input = ');
Max_Input = xlsread('data',1,'J2');
%Max_Output = input('\n最多几个输出端口：\nMax_Output = ');
Max_Output = xlsread('data',1,'L2');

%Max_InputNum = input('\n某个输入端口的，最多输入数量\nMax_InputNum = ');      %决定矩阵每一行元素数量
Max_InputNum = xlsread('data',1,'K2');
%Max_OutputNum = input('\n某个输出端口的，最多输出数量\nMax_OutputNum = ');     %决定矩阵每一行元素数量
Max_OutputNum = xlsread('data',1,'M2');

%********************************************************************************
%********************************************************************************
%********************************************************************************
%*****************************【生成Switch.xls】*************************************
global S_InSeq;               %设置输入序列
S_In = zeros(Max_Input,Max_InputNum,Num);  %行Max_Input代表最多几个输入端口，
                                           %列Max_InputNum代表某个输入端口最多的输入数量，
                                           %页Num代表交换机数量
S_Out = zeros(Max_Output,Max_OutputNum,Num);  %行代表最多几个输出端口，
                                              %列代表某个输出端口最多的输出数量，
                                              %页代表交换机数量
                                              
judge = exist('Switch.xls','file'); %判断是否存在旧的In.xls
if judge == 2   %如果目前存在已有In.xls文件
    select = input('是否需要删除已有Switch.xls并重新生成？ Y/N\n','s');
    if select == 'Y'||'y'    %需要删除
        delete('Switch.xls');   %删除Switch.xls以方便后续重新生成最新Switch.xls表格
        disp('已成功删除旧表！');
        %生成Switch.xls表格以键入数据
        disp('正在创建新表……');
        range_in = [num2str(2),':',num2str(Max_Input + 2 - 1)];
        range_out = [num2str(Max_Input + 4),':',num2str(Max_Output + Max_Input + 4 - 1)];
        for i = 1:Num
            eval(['xlswrite(''Switch'',S_In(:,:,',num2str(i),'),',num2str(i),',range_in);']);    %Sheet i就代表Si的输入数据.把【输入】的纯0矩阵写到Switch.xls中
            eval(['xlswrite(''Switch'',S_Out(:,:,',num2str(i),'),',num2str(i),',range_out);']);    %Sheet i就代表Si的输入数据.把【输出】的纯0矩阵写到Switch.xls中

            switch_in = {'这是输入矩阵'};
            switch_out = {'这是输出矩阵'};
            eval(['xlswrite(''Switch'',switch_in,',num2str(i),',''A',num2str(1),''');']);    %在第一行写“这是输入矩阵”
            eval(['xlswrite(''Switch'',switch_out,',num2str(i),',''A',num2str(Max_Input + 3),''');']);    %在输出矩阵的上一行写“这是输出矩阵”
            
            switch_loca = {'交换机位置'};
            eval(['xlswrite(''Switch'',switch_loca,',num2str(i),',''A',num2str(Max_Output + Max_Input + 8),''');']);    %把下面隔两行的位置处，设置输入交换机位置的格子(8 = + 4 - 1 + 5)
            eval(['xlswrite(''Switch'',',num2str(1),',',num2str(i),',''A',num2str(Max_Output + Max_Input + 9),''');']);  %预设每个位置都是1.(9 = + 4 - 1 + 6)
            
            %{
            旧版title
            title1 = {'Sheet'};
            title2 = {'代表S'};
            eval(['xlswrite(''Switch'',title1,',num2str(i),',''A',num2str(Max_Output + Max_Input + 13),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''B',num2str(Max_Output + Max_Input + 13),''');']);
            eval(['xlswrite(''Switch'',title2,',num2str(i),',''C',num2str(Max_Output + Max_Input + 13),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''D',num2str(Max_Output + Max_Input + 13),''');']);
            %}
            title = {['Sheet',num2str(i),'代表S',num2str(i)]};
            eval(['xlswrite(''Switch'',title,',num2str(i),',''A',num2str(Max_Input + 15),''');']);
        end
    end
else if judge == 0       %如果不存在Switch.xls
        %生成Switch.xls表格以键入数据    
        disp('正在生成Switch.xls表格……');
        range_in = [num2str(2),':',num2str(Max_Input + 2 - 1)];
        range_out = [num2str(Max_Input + 4),':',num2str(Max_Output + Max_Input + 4 - 1)];
        for i = 1:Num
            eval(['xlswrite(''Switch'',S_In(:,:,',num2str(i),'),',num2str(i),',range_in);']);    %Sheet i就代表Si的输入数据.把【输入】的纯0矩阵写到Switch.xls中
            eval(['xlswrite(''Switch'',S_Out(:,:,',num2str(i),'),',num2str(i),',range_out);']);    %Sheet i就代表Si的输入数据.把【输出】的纯0矩阵写到Switch.xls中

            switch_in = {'这是输入矩阵'};
            switch_out = {'这是输出矩阵'};
            eval(['xlswrite(''Switch'',switch_in,',num2str(i),',''A',num2str(1),''');']);    %在第一行写“这是输入矩阵”
            eval(['xlswrite(''Switch'',switch_out,',num2str(i),',''A',num2str(Max_Input + 3),''');']);    %在输出矩阵的上一行写“这是输出矩阵”
            
            switch_loca = {'交换机位置'};
            eval(['xlswrite(''Switch'',switch_loca,',num2str(i),',''A',num2str(Max_Output + Max_Input + 8),''');']);    %把下面隔两行的位置处，设置输入交换机位置的格子(8 = + 4 - 1 + 5)
            eval(['xlswrite(''Switch'',',num2str(1),',',num2str(i),',''A',num2str(Max_Output + Max_Input + 9),''');']);  %预设每个位置都是1.(9 = + 4 - 1 + 6)

            %{
            旧版title
            title1 = {'Sheet'};
            title2 = {'代表S'};
            eval(['xlswrite(''Switch'',title1,',num2str(i),',''A',num2str(Max_Input + 15),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''B',num2str(Max_Input + 15),''');']);
            eval(['xlswrite(''Switch'',title2,',num2str(i),',''C',num2str(Max_Input + 15),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''D',num2str(Max_Input + 15),''');']);
            %}
            title = {['Sheet',num2str(i),'代表S',num2str(i)]};
            eval(['xlswrite(''Switch'',title,',num2str(i),',''A',num2str(Max_Input + 15),''');']);
        end
    end
end

disp('请打开Switch.xls文件,并在Sheet_X中键入交换机S_X的数据，保存后请关闭表格，并敲击回车。');
disp('（若无需动作请键入回车）');
disp(' ');
pause;
disp('----------------------------');
disp('所有交换机输入数据，键入完成！');
disp('请稍后……正在处理中');
disp('----------------------------');
disp(' ');

%把Switch.xls表格内的【输入】内容读到S_In内。Sheet1为S1，Sheet2为S2
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
    %读取【输入】数据
    for j = 1:Max_Input
            range = [num2str(j + 1),':',num2str(j + 1)];    %读取第j行
            S_In(j,:,i) = xlsread('Switch',i,range);     %按一行一行读,得到真正S_In
    end
    %读取【输出】数据
    for j = 1:Max_Output
        range = [num2str(Max_Input + 3 + j),':',num2str(Max_Input + 3 + j)];    %读取第Max_Input + 3 + j行，其中3是因为“两个表格名字 + 一个空行”
        S_Out(j,:,i) = xlsread('Switch',i,range);     %按一行一行读,得到真正S_Out
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
    range_in = [num2str(2),':',num2str(Max_Input + 2 - 1)];
    range_out = [num2str(Max_Input + 4),':',num2str(Max_Output + Max_Input + 4 - 1)];
    
    temp_in(:,:) = xlsread('Switch',i,range_in);     %读取了S1、S2、S3的输入
    temp_out(:,:) = xlsread('Switch',i,range_out);
    seq = xlsread('Switch',i,['A',num2str(Max_Output + Max_Input + 9)]);     %seq代表这个交换机的位置（第几个）
    S_SwitchSeq(i) = seq;   %弄出S几是第几个交换机。
    S_InSeq(:,:,seq) = temp_in;
    S_OutSeq(:,:,seq) = temp_out;
end

for i = 1:Num
    n = 1;
    while ~isequal(S_In(:,:,n),S_InSeq(:,:,i))  %无论用S_In还是S_Out来比较都可以
        n = n + 1;
    end
    S_ReSwitchSeq(i) = n;   %弄出第几个交换机是S几。
end

%********************************************************************************
%********************************************************************************
%********************************************************************************
%********************************************************************************

disp('***********************【网络输入分布矩阵图】***********************');
S_InSeq
disp(' ');
disp('***********************【网络输出分布矩阵图】***********************');
S_OutSeq
clear judge select switch_loca title1 title2 range temp seq i j n;


%% **************【整理S_Out与S_OutSeq】***********************
%Max_Output = input('\n最多几个输出端口：\nMax_Output = ');
Max_Output = xlsread('data',1,'L2');

%Max_OutputNum = input('\n某个输出端口的，最多输出数量\nMax_OutputNum = ');     %决定矩阵每一行元素数量
Max_OutputNum = xlsread('data',1,'M2');
disp(['下列矩阵元素个数为',num2str(Max_OutputNum),'个']);

global S_OutSeq;               %设置输出序列
S_Out = zeros(Max_Output,Max_OutputNum,Num);  %行代表最多几个输出端口，
                                              %列代表某个输出端口最多的输出数量，
                                              %页代表交换机数量
   for i = 1:Num        %为了得到每一个交换机的输出序列
       disp('*******************************************');
       for j = 1:Max_Output
           disp(' ');
           disp(['请键入S',num2str(i),'交换机的输出数据，']);
           disp('数据键入格式为行向量，保证“元素个数 = 最多输出数量”，如：[a b c]。');
           S_Out(j,:,i) = input(' ');
           disp(' ');
       end
   end
clear i j;
%实现：第1页是S3，第2页是S1，第3页是S2。
for i = 1:Num
    disp(['这是第',num2str(i),'个交换机，请确认这是S几，并将这个数字代替x：S_Out(:,:,x)']);
    x=input('');
    temp = S_Out(:,:,x);
    S_OutSeq(:,:,i) = temp;       %S_OutSeq是具有拓扑顺序的交换机结构，其第一页就代表第一个交换机
    S_OutSeqNo(i) = x;            %S_OutSeqNo代表第几个交换机是S几，可以直接得到序号
    disp(' ');
end
clear x temp i;
%不清除S_Out是因为要留下直接用S几来确定输入输出参数