%����Ҫ�ϲ�S_In��S_Out

%% S_In
%ÿһ�д���һ������˿ڣ�ÿһ��û��ʵ�ʺ��壨һ���е�ĳһ�д�������˿��е�һ��֡��������ÿһҳ����һ��������

%Num = input('\n������������\nNum = ');
Num = xlsread('data',1,'I2');

%Max_Input = input('\n��༸������˿ڣ�\nMax_Input = ');
Max_Input = xlsread('data',1,'J2');
%Max_Output = input('\n��༸������˿ڣ�\nMax_Output = ');
Max_Output = xlsread('data',1,'L2');

%Max_InputNum = input('\nĳ������˿ڵģ������������\nMax_InputNum = ');      %��������ÿһ��Ԫ������
Max_InputNum = xlsread('data',1,'K2');
%Max_OutputNum = input('\nĳ������˿ڵģ�����������\nMax_OutputNum = ');     %��������ÿһ��Ԫ������
Max_OutputNum = xlsread('data',1,'M2');

%********************************************************************************
%********************************************************************************
%********************************************************************************
%*****************************������Switch.xls��*************************************
global S_InSeq;               %������������
S_In = zeros(Max_Input,Max_InputNum,Num);  %��Max_Input������༸������˿ڣ�
                                           %��Max_InputNum����ĳ������˿���������������
                                           %ҳNum������������
S_Out = zeros(Max_Output,Max_OutputNum,Num);  %�д�����༸������˿ڣ�
                                              %�д���ĳ������˿��������������
                                              %ҳ������������
                                              
judge = exist('Switch.xls','file'); %�ж��Ƿ���ھɵ�In.xls
if judge == 2   %���Ŀǰ��������In.xls�ļ�
    select = input('�Ƿ���Ҫɾ������Switch.xls���������ɣ� Y/N\n','s');
    if select == 'Y'||'y'    %��Ҫɾ��
        delete('Switch.xls');   %ɾ��Switch.xls�Է������������������Switch.xls���
        disp('�ѳɹ�ɾ���ɱ�');
        %����Switch.xls����Լ�������
        disp('���ڴ����±���');
        range_in = [num2str(2),':',num2str(Max_Input + 2 - 1)];
        range_out = [num2str(Max_Input + 4),':',num2str(Max_Output + Max_Input + 4 - 1)];
        for i = 1:Num
            eval(['xlswrite(''Switch'',S_In(:,:,',num2str(i),'),',num2str(i),',range_in);']);    %Sheet i�ʹ���Si����������.�ѡ����롿�Ĵ�0����д��Switch.xls��
            eval(['xlswrite(''Switch'',S_Out(:,:,',num2str(i),'),',num2str(i),',range_out);']);    %Sheet i�ʹ���Si����������.�ѡ�������Ĵ�0����д��Switch.xls��

            switch_in = {'�����������'};
            switch_out = {'�����������'};
            eval(['xlswrite(''Switch'',switch_in,',num2str(i),',''A',num2str(1),''');']);    %�ڵ�һ��д�������������
            eval(['xlswrite(''Switch'',switch_out,',num2str(i),',''A',num2str(Max_Input + 3),''');']);    %������������һ��д�������������
            
            switch_loca = {'������λ��'};
            eval(['xlswrite(''Switch'',switch_loca,',num2str(i),',''A',num2str(Max_Output + Max_Input + 8),''');']);    %����������е�λ�ô����������뽻����λ�õĸ���(8 = + 4 - 1 + 5)
            eval(['xlswrite(''Switch'',',num2str(1),',',num2str(i),',''A',num2str(Max_Output + Max_Input + 9),''');']);  %Ԥ��ÿ��λ�ö���1.(9 = + 4 - 1 + 6)
            
            %{
            �ɰ�title
            title1 = {'Sheet'};
            title2 = {'����S'};
            eval(['xlswrite(''Switch'',title1,',num2str(i),',''A',num2str(Max_Output + Max_Input + 13),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''B',num2str(Max_Output + Max_Input + 13),''');']);
            eval(['xlswrite(''Switch'',title2,',num2str(i),',''C',num2str(Max_Output + Max_Input + 13),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''D',num2str(Max_Output + Max_Input + 13),''');']);
            %}
            title = {['Sheet',num2str(i),'����S',num2str(i)]};
            eval(['xlswrite(''Switch'',title,',num2str(i),',''A',num2str(Max_Input + 15),''');']);
        end
    end
else if judge == 0       %���������Switch.xls
        %����Switch.xls����Լ�������    
        disp('��������Switch.xls��񡭡�');
        range_in = [num2str(2),':',num2str(Max_Input + 2 - 1)];
        range_out = [num2str(Max_Input + 4),':',num2str(Max_Output + Max_Input + 4 - 1)];
        for i = 1:Num
            eval(['xlswrite(''Switch'',S_In(:,:,',num2str(i),'),',num2str(i),',range_in);']);    %Sheet i�ʹ���Si����������.�ѡ����롿�Ĵ�0����д��Switch.xls��
            eval(['xlswrite(''Switch'',S_Out(:,:,',num2str(i),'),',num2str(i),',range_out);']);    %Sheet i�ʹ���Si����������.�ѡ�������Ĵ�0����д��Switch.xls��

            switch_in = {'�����������'};
            switch_out = {'�����������'};
            eval(['xlswrite(''Switch'',switch_in,',num2str(i),',''A',num2str(1),''');']);    %�ڵ�һ��д�������������
            eval(['xlswrite(''Switch'',switch_out,',num2str(i),',''A',num2str(Max_Input + 3),''');']);    %������������һ��д�������������
            
            switch_loca = {'������λ��'};
            eval(['xlswrite(''Switch'',switch_loca,',num2str(i),',''A',num2str(Max_Output + Max_Input + 8),''');']);    %����������е�λ�ô����������뽻����λ�õĸ���(8 = + 4 - 1 + 5)
            eval(['xlswrite(''Switch'',',num2str(1),',',num2str(i),',''A',num2str(Max_Output + Max_Input + 9),''');']);  %Ԥ��ÿ��λ�ö���1.(9 = + 4 - 1 + 6)

            %{
            �ɰ�title
            title1 = {'Sheet'};
            title2 = {'����S'};
            eval(['xlswrite(''Switch'',title1,',num2str(i),',''A',num2str(Max_Input + 15),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''B',num2str(Max_Input + 15),''');']);
            eval(['xlswrite(''Switch'',title2,',num2str(i),',''C',num2str(Max_Input + 15),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''D',num2str(Max_Input + 15),''');']);
            %}
            title = {['Sheet',num2str(i),'����S',num2str(i)]};
            eval(['xlswrite(''Switch'',title,',num2str(i),',''A',num2str(Max_Input + 15),''');']);
        end
    end
end

disp('���Switch.xls�ļ�,����Sheet_X�м��뽻����S_X�����ݣ��������رձ�񣬲��û��س���');
disp('�������趯�������س���');
disp(' ');
pause;
disp('----------------------------');
disp('���н������������ݣ�������ɣ�');
disp('���Ժ󡭡����ڴ�����');
disp('----------------------------');
disp(' ');

%��Switch.xls����ڵġ����롿���ݶ���S_In�ڡ�Sheet1ΪS1��Sheet2ΪS2
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
    %��ȡ�����롿����
    for j = 1:Max_Input
            range = [num2str(j + 1),':',num2str(j + 1)];    %��ȡ��j��
            S_In(j,:,i) = xlsread('Switch',i,range);     %��һ��һ�ж�,�õ�����S_In
    end
    %��ȡ�����������
    for j = 1:Max_Output
        range = [num2str(Max_Input + 3 + j),':',num2str(Max_Input + 3 + j)];    %��ȡ��Max_Input + 3 + j�У�����3����Ϊ������������� + һ�����С�
        S_Out(j,:,i) = xlsread('Switch',i,range);     %��һ��һ�ж�,�õ�����S_Out
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
    range_in = [num2str(2),':',num2str(Max_Input + 2 - 1)];
    range_out = [num2str(Max_Input + 4),':',num2str(Max_Output + Max_Input + 4 - 1)];
    
    temp_in(:,:) = xlsread('Switch',i,range_in);     %��ȡ��S1��S2��S3������
    temp_out(:,:) = xlsread('Switch',i,range_out);
    seq = xlsread('Switch',i,['A',num2str(Max_Output + Max_Input + 9)]);     %seq���������������λ�ã��ڼ�����
    S_SwitchSeq(i) = seq;   %Ū��S���ǵڼ�����������
    S_InSeq(:,:,seq) = temp_in;
    S_OutSeq(:,:,seq) = temp_out;
end

for i = 1:Num
    n = 1;
    while ~isequal(S_In(:,:,n),S_InSeq(:,:,i))  %������S_In����S_Out���Ƚ϶�����
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
disp(' ');
disp('***********************����������ֲ�����ͼ��***********************');
S_OutSeq
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