clear all;
clc;

%% *************************[The FORMULA]*****************************
%{
  W_lasti_i_t = D_hpX_i_t + D_Xmax_i_t + D_lpX_i + sigma_max_Cj + ...
    (num(Pi)-1)*sl - positive(sigma_delta_h_i_t - t) - Ci*(1+max(alfa_neg_X / alfa_pos_X))

  Ri = max(W_lasti_i_t - t + Ci)
%}


%% ************************[Show the Network Diagram]*******************
%map=imread('topology_eg.jpg');
%imshow(map);


%% *****************[Define the Transmission Parameters according to Model]***********************
%{
disp('Attention! ');
disp('The data input below should be matrix or vector! ');
%Tao = input('\nPlease enter the number of flows: \n�� = ');
%}
Tao = xlsread('data',1,'A2');
global I;
I = xlsread('data',1,'C2');  % I is the subscript of the researched flow, which is the location of it in the flows.
%disp(' ');
%disp('***********************************');
global C;
for i = 1:Tao
           %disp(['Please enter C, which is the transmission time of each frame',num2str(i),'/total number is',num2str(Tao),'(ms)']);
           %C(i) = input('C = ');
           %eval(['C(i) = xlsread(','''data''',num2str(1),'''F',num2str(i+1),''')']);
           eval(['C(i) = xlsread(''data'',1,','''F',num2str(i+1),''');']);
       end
       C
%{
       for i = 1:Tao
           disp(['Please enter C, which is the transmission time of each frame',num2str(i),'/total number is',num2str(Tao),'(ms)']);
           C(i) = input('C = ');
           disp(' ');
       end

       %C=[40 40 40 40 40 120];                 %transmission time(meaning the length of frame)
       %C = input('\nPlease enter C, which is the transmission time of each frame(ms)��\nC = ');
disp('***********************************');
%}

global T;
for i = 1:Tao
           %disp(['Please enter T, which is the period of sending of each frame from ECU.',num2str(i),'/the total number is',num2str(Tao),'(ms)']);
           %T(i) = input('T = ');
           eval(['T(i) = xlsread(''data'',1,','''G',num2str(i+1),''');']);
       end
       T
%{
       for i = 1:Tao
           disp(['Please enter T, which is the period of sending of each frame from ECU.',num2str(i),'/the total number is',num2str(Tao),'(ms)']);
           T(i) = input('T = ');
           disp(' ');
       end
       %T=[2000 2000 2000 160 2000 2000];       %sending period
       %T = input('\nPlease enter T, which is the period of sending of each frame from ECU.(ms)��\nT = ');
disp('***********************************');
%}
global Pr;
for i = 1:Tao
           %disp(['Please enter Pr, which is the priority of each flow, from��1 to ��n',num2str(i),'/the totla number is',num2str(Tao),'(1/2/3)']);
           %Pr_temp = input('Pr = ');
           eval(['Pr_temp = xlsread(''data'',1,','''H',num2str(i+1),''');']);
           while (Pr_temp ~= 1)&&(Pr_temp ~= 2)&&(Pr_temp ~= 3)         %check if every input is the element of {1,2,3}
                 disp(['Please enter the data of No.',num2str(i),'! The priority only can be the element of {1, 2, 3}! ']);
                 Pr_temp = input('Pr = ');
           end
           Pr(i) = Pr_temp;
       end
       Pr 

if ismember(3,Pr)       % if there is any Class A flow
    %{
   "the code below is replaced by function find"
    A = 1;
    for i=1:length(Pr)
        if Pr(i) == 3       %the priority of Class A is the highest one, which equals 3
            j_a(A) = i;     %pick up the Class A from Pr
        end
            A = A + 1;
    end
    A = 1;
    for i=1:length(j_a)
        if j_a(i) ~= 0
            temp = j_a(i);
            J_A(A) = temp;      %pick up the non-zero element from j_a, J_A is the set of all Class A flows.
            A = A + 1;
        end
    end
    
    A = 1;
    for i = 1:length(j_a)
        if (j_a(i) ~= 0) && (j_a(i) ~= I)       %pick up the non-I intefere flows from J_A
            temp = j_a(i);
            j_Reala(A) = temp;
        end
        A = A + 1;
    end
    A = 1;
    for i=1:length(j_Reala)
        if j_Reala(i) ~= 0
            temp = j_Reala(i);
            J_RealA(A) = temp;      %pick up the non-zero element from j_Reala, J_RealA is the flows from Class A apart from the researched flow.  //[4]
            A = A + 1;              %is this step is at the out of "if", then it cannot remove the 0 element.
        end
    end
    %}
    J_A = find(Pr == 3);
    J_A_temp = find(Pr == 3);
    %if ismember(I,J_A)      %if I belongs to Class A
    J_A_temp(find(J_A_temp == I)) = [];   %delete the data at the I's location
        %location_I = find(J_A == I);
        %J_A(location_I) = [];
    J_RealA = J_A_temp;      %J_RealA must exist
    %end
    clear location_I J_A_temp;
end

if ismember(2,Pr)       %if there is any Class B flow
    %{
    "the code below is replaced by function find"
    B = 1;
    for i=1:length(Pr)
        if Pr(i) == 2       %the priority of Class B equals 2
            j_b(B) = i;     %pick up the Class B from Pr
        end
        B = B + 1;
    end
    B = 1;
    for i=1:length(j_b)
        if j_b(i) ~= 0
            temp = j_b(i);
            J_B(B) = temp;      %pick up the non-zero element from j_b, J_B is the set of all Class A flows.  //[1 2 3]
            B = B + 1;
        end
    end
    
    B = 1;
    for i = 1:length(j_b)
        if (j_b(i) ~= 0) && (j_b(i) ~= I)
            temp = j_b(i);
            j_Realb(B) = temp;  %remove the researched flow I which has the same Class with the interfere flows. 
        end
        B = B + 1;
    end
    B = 1;
    for i=1:length(j_Realb)
        if j_Realb(i) ~= 0
            temp = j_Realb(i);
            J_RealB(B) = temp;      %pick up the non-zero element from j_Realb, J_RealB is the flows from Class B apart from the researched flow.  //[2 3]
            B = B + 1;              %is this step is at the out of "if", then it cannot remove the 0 element.
        end
    end
    %}
    J_B = find(Pr == 2);
    J_B_temp = find(Pr == 2);
    %if ismember(I,J_B)      %if I belongs to Class B
    J_B_temp(find(J_B_temp == I)) = [];
        %location_I = find(J_B == I);
        %J_B(location_I) = [];
    J_RealB = J_B_temp;
    %end
    clear location_I J_B_temp;
end
       
disp('***************************');
R = xlsread('data',1,'B2');    % Mbit/s
%R = input('\ntransmission rate:\nR = ');
Alfa = R * [0.5 0.5 0.25 0.75];
global alfa_neg_A;
global alfa_pos_A;
global alfa_neg_B;
global alfa_pos_B;
       alfa_pos_A = Alfa(1);
       alfa_neg_A = -Alfa(2);
       alfa_pos_B = Alfa(3);
       alfa_neg_B = -Alfa(4);

      
%I = input('\nPlease enter the subscript of researched flow\nI = ');
global sl;
sl = xlsread('data',1,'D2');       
%sl = input('\nPlease enter switch latency\nsl = ');
global t;
t = xlsread('data',1,'E2');
%t=input('\nPlease enter t, which is the generation time of I. \n');

%I cannot compete with itself. 
S_first_ii_max_i = 0;
S_last_ii_min_i = 0;
A_ii = 0;
%*************************************************************


%% ************[Get the S_In, S_Out, S_InSeq and S_OutSeq]*****************
%each row represents an input port.
%each colomn does not have a real meaning.(any colomn of one of rows is a frame or flow from input port)
%each page means a switch machine. 
%{
Num = input('\nthe number of switch machine: \nNum = ');
Max_Input = input('\nthe maximum number of input port: \nMax_Input = ');
Max_InputNum = input('\nthe maximum number of input flows: \nMax_InputNum = ');      %decide the number of colomns
%}

%Num = input('\nthe number of switch machine: \nNum = ');
Num = xlsread('data',1,'I2');

%Max_Input = input('\nthe maximum number of input port: \nMax_Input = ');
Max_Input = xlsread('data',1,'J2');
%Max_Output = input('\nthe maximum number of output port: \nMax_Output = ');
Max_Output = xlsread('data',1,'L2');

%Max_InputNum = input('\nthe maximum number of input flows: \nMax_InputNum = ');      %decide the number of rows
Max_InputNum = xlsread('data',1,'K2');
%Max_OutputNum = input('\nthe maximum number of output flows: \nMax_OutputNum = ');     %��������ÿһ��Ԫ������
Max_OutputNum = xlsread('data',1,'M2');

%********************************************************************************
%********************************************************************************
%********************************************************************************
%*****************************[Generate Switch_test.xls]*************************************
global S_InSeq;               %Set the input series
global S_OutSeq;
S_In = zeros(Max_Input,Max_InputNum,Num);  %the Row, Max_Input, means the maximum number of input port; 
                                           %the colomn, Max_InputNum, means the maximum number of input flows; 
                                           %the page, Num, means the number of switch machine. 
S_Out = zeros(Max_Output,Max_OutputNum,Num);  %the Row, Max_Input, means the maximum number of input port; 
                                              %the colomn, Max_InputNum, means the maximum number of input flows; 
                                              %the page, Num, means the number of switch machine. 

%decide if it is necessary to delete the old Switch_test.xls, and create a new one.
judge = exist('Switch_test.xls','file');
if judge == 2   %if there does exist Switch_test.xls
    select = input('Do you need to delete the Switch_test.xls and recreate one? Y/N\n','s');
    if select == 'Y'&&'y'
        delete('Switch_test.xls');   %delete Switch.xls to generate the latest Switch.xls
        disp('Delete the old sheet successefully! ');
        %generate Switch.xls to enter the data
        disp('*************************');
        disp('Creating new sheet...');
        disp('*************************');
        range_in = [num2str(2),':',num2str(Max_Input + 2 - 1)];
        range_out = [num2str(Max_Input + 4),':',num2str(Max_Output + Max_Input + 4 - 1)];
        for i = 1:Num
            eval(['xlswrite(''Switch_test'',S_In(:,:,',num2str(i),'),',num2str(i),',range_in);']);    %Sheet_i represents the data of S_i, and write the 0 matrix of input into Switch_test.xls
            eval(['xlswrite(''Switch_test'',S_Out(:,:,',num2str(i),'),',num2str(i),',range_out);']);    %Sheet_i represents the data of S_i, and write the 0 matrix of output into Switch_test.xls

            switch_in = {'Input'};
            switch_out = {'Output'};
            eval(['xlswrite(''Switch_test'',switch_in,',num2str(i),',''A',num2str(1),''');']);    %write "Input" at row 1.
            eval(['xlswrite(''Switch_test'',switch_out,',num2str(i),',''A',num2str(Max_Input + 3),''');']);    %write "Output" at the row on the output matrix
            
            switch_loca = {'Location of Switch'};
            eval(['xlswrite(''Switch_test'',switch_loca,',num2str(i),',''A',num2str(Max_Output + Max_Input + 8),''');']);    %the interval is 2 rows(8 = + 4 - 1 + 5)
            eval(['xlswrite(''Switch_test'',',num2str(1),',',num2str(i),',''A',num2str(Max_Output + Max_Input + 9),''');']);  %pre-set each value is 1.(9 = + 4 - 1 + 6)
            
            %{
            the old version of title
            title1 = {'Sheet'};
            title2 = {'represents S'};
            eval(['xlswrite(''Switch'',title1,',num2str(i),',''A',num2str(Max_Output + Max_Input + 13),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''B',num2str(Max_Output + Max_Input + 13),''');']);
            eval(['xlswrite(''Switch'',title2,',num2str(i),',''C',num2str(Max_Output + Max_Input + 13),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''D',num2str(Max_Output + Max_Input + 13),''');']);
            %}
            title = {['Sheet',num2str(i),'represents S',num2str(i)]};
            eval(['xlswrite(''Switch_test'',title,',num2str(i),',''A',num2str(Max_Input + 15),''');']);
        end
    end
else if judge == 0       %if there is not Switch_test.xls
        %generate Switch.xls
        disp('********************************');
        disp('Creating Switch_test.xls ...');
        disp('********************************');
        range_in = [num2str(2),':',num2str(Max_Input + 2 - 1)];
        range_out = [num2str(Max_Input + 4),':',num2str(Max_Output + Max_Input + 4 - 1)];
        for i = 1:Num
            eval(['xlswrite(''Switch_test'',S_In(:,:,',num2str(i),'),',num2str(i),',range_in);']);    %Sheet_i represents the data of S_i, and write the 0 matrix of input into Switch_test.xls
            eval(['xlswrite(''Switch_test'',S_Out(:,:,',num2str(i),'),',num2str(i),',range_out);']);    %Sheet_i represents the data of S_i, and write the 0 matrix of output into Switch_test.xls

            switch_in = {'Input'};
            switch_out = {'Output'};
            eval(['xlswrite(''Switch_test'',switch_in,',num2str(i),',''A',num2str(1),''');']);    %write "Input" at row 1.
            eval(['xlswrite(''Switch_test'',switch_out,',num2str(i),',''A',num2str(Max_Input + 3),''');']);    %write "Output" at the row on the output matrix
            
            switch_loca = {'Switch Location'};
            eval(['xlswrite(''Switch_test'',switch_loca,',num2str(i),',''A',num2str(Max_Output + Max_Input + 8),''');']);    %the interval is 2 rows(8 = + 4 - 1 + 5)
            eval(['xlswrite(''Switch_test'',',num2str(1),',',num2str(i),',''A',num2str(Max_Output + Max_Input + 9),''');']);  %pre-set each value is 1.(9 = + 4 - 1 + 6)

            %{
            the old version of title
            title1 = {'Sheet'};
            title2 = {'represents S'};
            eval(['xlswrite(''Switch'',title1,',num2str(i),',''A',num2str(Max_Input + 15),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''B',num2str(Max_Input + 15),''');']);
            eval(['xlswrite(''Switch'',title2,',num2str(i),',''C',num2str(Max_Input + 15),''');']);
            eval(['xlswrite(''Switch'',',num2str(i),',',num2str(i),',''D',num2str(Max_Input + 15),''');']);
            %}
            title = {['Sheet',num2str(i),'represents S',num2str(i)]};
            eval(['xlswrite(''Switch_test'',title,',num2str(i),',''A',num2str(Max_Input + 15),''');']);
        end
    end
end
disp('If you need to modify some parameters, open Switch_test.xls.');
disp('Enter the data of switch S_X into each Sheet_X. ');
disp('Save and close the sheet, and click the Enter. ');
disp('--------------------------------------------------------------------');
disp('(If you do not need to modify any parameter, just click the Enter. )');
disp(' ');
pause;
disp('****************************************************');
disp('Processing all of the data in the DATA SHEET...');
disp('Please wait...');
disp('****************************************************');
disp(' ');

%read the input data from Switch_test.xls to S_In. Sheet1 is S1��Sheet2 is S2
%the old versioin: 
%{
    "the code below is the old version"
       for i = 1:Num
           disp('********************');
           for j = 1:Max_Input
               disp(['Please enter the data of S',num2str(i),'and the port ',num2str(j),'. ']);
               disp(['data format is line vector, and keep that the number of elements equals the Max_InputNum��eg. [a b c]']);
               S_In(j,:,i) = input(' ');
               disp(' ');
           end
       end
%}  
%the new version: 
for i = 1:Num
    %read the INPUT data
    for j = 1:Max_Input
            range = [num2str(j + 1),':',num2str(j + 1)];    %read the row j
            S_In(j,:,i) = xlsread('Switch_test',i,range);     %read row by row, and get the real S_In
    end
    %read the OUTPUT data
    for j = 1:Max_Output
        range = [num2str(Max_Input + 3 + j),':',num2str(Max_Input + 3 + j)];    %read the row of No.(Max_Input + 3 + j), in which 3 means that: two names of sheet and a blank row
        S_Out(j,:,i) = xlsread('Switch_test',i,range);     %read row by row, and get the real S_Out
    end
    
end


%S1=[1 0 0;2 3 6;0 0 0];
%S2=[1 2 3;4 0 0;5 0 0];
%S3=[2 0 0;3 0 0;6 0 0];

%Intention: the 1st page is S3, the 2nd page is S1, the 3rd page is S2. And get S_InSeq. 
%the old version: 
%{
OLD:
for i = 1:Num
    disp(['This is the No.',num2str(i),', and confirm this is S_x, and replace the x with the specific number: S_In(:,:,x)']);
    x=input('');
    temp = S_In(:,:,x);
    S_InSeq(:,:,i) = temp;       %S_InSeq is a structure with topology of the network, its first page represents the first switch machine.
    disp(' ');
end
%}
%the new version: 
for i = 1:Num
    range_in = [num2str(2),':',num2str(Max_Input + 2 - 1)];
    range_out = [num2str(Max_Input + 4),':',num2str(Max_Output + Max_Input + 4 - 1)];
    
    temp_in(:,:) = xlsread('Switch_test',i,range_in);     %read the input of S1, S2 and S3
    temp_out(:,:) = xlsread('Switch_test',i,range_out);
    seq = xlsread('Switch_test',i,['A',num2str(Max_Output + Max_Input + 9)]);     %seq represents the location of this swtich
    S_SwitchSeq(i) = seq;   %get clear the location of S_x
    S_InSeq(:,:,seq) = temp_in;
    S_OutSeq(:,:,seq) = temp_out;
end

for i = 1:Num
    n = 1;
    while ~isequal(S_In(:,:,n),S_InSeq(:,:,i))  %no matter compared with S_In or S_Out
        n = n + 1;
    end
    S_ReSwitchSeq(i) = n;   %get clear that any switch is S_which. 
    S_OutSeqNo(i) = n;      %the real meaning is the same as S_ReSwitchSeq, 
                            %but do not need to modify other parameters. 
end

%********************************************************************************
%********************************************************************************
%********************************************************************************
%********************************************************************************

disp('***********************[the INPUT matrix]***********************');
S_InSeq
disp(' ');
disp('***********************[the OUTPUT matrix]***********************');
S_OutSeq
clear judge select switch_loca title1 title2 range temp seq i j n Ifopen;


%% ***************[Search for Pi, the Path of I]************************
%the path matrix,
%first row means the S_X,
%each colomn means an output port at a switch
p = 0;
for i = 1:Num
    for j = 1:Max_Output
        m = 0;
        for k = 1:Max_OutputNum
            if S_OutSeq(j,k,i) == I     %I exists, meaning this output port belongs to Pi
                m=m+1;
                p=p+1;
            end
            if m == 1
                Pi(:,p) = [S_OutSeqNo(i),j];        %"Pi" is path matrix(the first row is switch number)
            end
        end
    end
end
clear i j m p;


%% **********************[Find source_i/j]************************************
%Find every flow's birth place.
%See the switch as origin. Eg. flow 2 was generated at S3, 
%and S3 is the 1st switch, so the birth place of flow 2 is 1.
source = zeros(Tao,Num);    %Total number of flows is Tao, and the maximum number of switch with input is Num,
                            %which is avoid that last switch only transmit one data to next switch
count_source = 1;
for i = 1:Num
    for j = 1:Max_Input
        n = 0;
        for k = 1:Max_InputNum
            if S_In(j,k,i) ~= 0
                n = n + 1;
            end
        end
        if n == 1    %this input port only has one input
            loca_source = find(S_In(j,:,i) ~= 0);    %find the location of single input at every input port(avoid the random input of human action)
            source_flow = S_In(j,loca_source,i);     %the subscript of single input at input port
            count_source = 1;
            while ~isequal(S_In(:,:,i),S_InSeq(:,:,count_source))     
                count_source = count_source + 1;      
            end
            source(source_flow,i) = count_source;     %"source" matrix means: 
                                                      %the element means the Nth switch
                                                      %the number of colomn means S_X
                                                      %the number of row means subscript of flow
        end
    end
end
for i = 1:Tao
    temp = find(source(i,:) ~= 0);
    Flow_Source(i) = source(i,temp);    %Flow_Source represents the start location of each flow. 
end
clear source source_flow loca_source count_source i temp;


%% **********************[Find first_ij]*************************************
%scan all of the output belongs to Pi, 
%Sequent the switches own the i and j in S_InSeq,
%get the first one out, and this is first_ij,
%the last one is last_ij. 
%Chinese Intention: �Ѿ���i��j�Ľ�������S_InSeq���Ÿ��򣬰ѵ�һ���ó�������first_ij�����һ������last_ij
disp(' ');
disp('********************************************');
disp('[Please note]: first_ij and last_ij: ');
%Find all Class B flows
if ismember(2,Pr)
    judge_first = 0;                %the parameter of judging if it is first_ij
    for J = 1:length(J_B)
        [hang lie] = size(Pi);
        for i = 1:lie        %from the first Switch of Pi to the last Switch of Pi
            %find out which is the switch when i=1
            Switch_No = Pi(1,i);    %the [i]th switch at Pi is S[Switch_No], which means that the [Switch_No]th page in S_Out
            for j = 1:Max_Output
                m = 0;
                %{
                Because of lacking the approaches to move out the influence of first_11[I = J_B(J)]
                
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
                if m == 2                % m==2 means that the output has I and J_B(J). ˵���������м���I��Ҳ��J_B(J),��ѭ���ֵ��ĸ�����
                    %cmpt(n) = i;         %(cmpt means compete, represents the [i]th switch)
                    %cmpt_out(n) = j;     %(cmpt_out represents the [j]th output port with compete at one of switches)
                    %n = n + 1;
                    m = 0;
                    judge_first = judge_first + 1;      %judge_first pluses itself with 1
                end
            end
            if judge_first == 1;                        %when accumulate the first time, meaning the first output has I and J_B, is first_ij
                first_ij = Switch_No;                   %represents the number of S
                for l = 1:Max_Output
                    if ismember(I,S_Out(l,:,first_ij))
                        first_ij_port = l;
                    end
                end
                disp(['first_',num2str(I),num2str(J_B(J)),' = S',num2str(first_ij),num2str(first_ij_port),', and the switch is S',num2str(first_ij)]);
                %store first_ij with structure Cmpt.first
                %query first_ij
                %without considering I = J_B(J)
                Cmpt.first(I,J_B(J),1) = first_ij;  %row means I, colomn means j, the first page means the number of switches, the second page means port number.
                Cmpt.first(I,J_B(J),2) = first_ij_port;
                break
            end
        end
        judge_first = 0;        %when finish passing through Pi, make judge_first zero. and change to the next priority of interfere. 
    end
    clear i j k m l first_ij_port;
end

%Find all the Class A flows
if ismember(3,Pr)
    judge_first = 0;                %the parameter of judging if it is first_ij
    for J = 1:length(J_A)
        [hang lie] = size(Pi);
        for i = 1:lie        %from the first Switch of Pi to the last Switch of Pi
            %find out which is the switch when i=1
            Switch_No = Pi(1,i);    %the [i]th switch at Pi is S[Switch_No], which means that the [Switch_No]th page in S_Out
            for j = 1:Max_Output
                m = 0;
                if ismember(J_A(J),S_Out(j,:,Switch_No))
                    for k = 1:Max_OutputNum
                        if (S_Out(j,k,Switch_No) == I) || (S_Out(j,k,Switch_No) == J_A(J))
                            m=m+1;
                        end
                    end
                end
                if m == 2                % m==2 means that the output has I and J_A(J). ˵���������м���I��Ҳ��J_A(J),��ѭ���ֵ��ĸ�����
                    %cmpt(n) = i;         %(cmpt means compete, represents the [i]th switch)
                    %cmpt_out(n) = j;     %(cmpt_out represents the [j]th output port with compete at one of switches)
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
                disp(['first_',num2str(I),num2str(J_A(J)),' = S',num2str(first_ij),num2str(first_ij_port),', and the switch is S',num2str(first_ij)]);
                %store first_ij with structure Cmpt.first
                %query first_ij
                %without considering I = J_A(J)
                Cmpt.first(I,J_A(J),1) = first_ij;  %row means I, colomn means j, the first page means the number of switches, the second page means port number.
                Cmpt.first(I,J_A(J),2) = first_ij_port;
                break
            end
        end
        judge_first = 0;
    end
    clear i j k m J judge_first first_ij Switch_No l first_ij_port hang lie;
end


%% *******************[Find last_ij]**************************
%{
������ͨ����ĳ����ͬʱ����I&J���ɴ˽������״̬��
�����߼���ͬʱ���ڣ����ж�Ϊ�棬ͨ����
��������ʧ������һ��������Ϊlast��
Ȼ���ٲ�һ��������һ�����������ĸ��˿ڡ�
���Ǽ���״̬�Ľ�����˳�����ʵ�֣�

New idea: find the first one with a reverse path
%}
reverse_Pi = fliplr(Pi);
%Find all the Class B flows
if ismember(2,Pr)
    judge_last = 0;                 %decide whether this is last_ij
    for J = 1:length(J_B)
        for i = 1:length(reverse_Pi)
            Switch_No = reverse_Pi(1,i);    %the [i]th switch at reverse_Pi is S[Switch_No], which means that the [Switch_No]th page in S_Out
            for j = 1:Max_Output    %port number
                m = 0;
                if ismember(J_B(J),S_Out(j,:,Switch_No))
                    for k = 1:Max_OutputNum
                        if (S_Out(j,k,Switch_No) == I) || (S_Out(j,k,Switch_No) == J_B(J))
                            m=m+1;
                        end
                    end
                end
                if m == 2      % m==2 means that the output has I and J_B(J). ˵���������м���I��Ҳ��J_B(J),��ѭ���ֵ��ĸ�����
                    %{
                    %[start monitoring last_ij, if these two flows keep existing, then it is "true", pass. ]
                    %[if they do not exist at the same time, then the last one is the "last", and test one more time for port]
                    whatchdog = 1;      %start monitoring
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
                disp(['last_',num2str(I),num2str(J_B(J)),' = S',num2str(last_ij),num2str(last_ij_port),', and the switch is S',num2str(last_ij)]);
                %store first_ij with structure Cmpt.first
                %query first_ij
                %without considering I = J_B(J)
                Cmpt.last(I,J_B(J),1) = last_ij;
                Cmpt.last(I,J_B(J),2) = last_ij_port;
                break
            end
        end
        judge_last = 0;
    end
    clear J i j m k l judge_last last_ij_port;
end

%Find all the Class A flows
if ismember(3,Pr)
    judge_last = 0;                         %decide whether this is last_ij
    for J = 1:length(J_A)
        [hang lie] = size(reverse_Pi);
        for i = 1:lie
            Switch_No = reverse_Pi(1,i);    %the [i]th switch at reverse_Pi is S[Switch_No], which means that the [Switch_No]th page in S_Out
            for j = 1:Max_Output    %port number
                m = 0;
                if ismember(J_A(J),S_Out(j,:,Switch_No))
                    for k = 1:Max_OutputNum
                        if (S_Out(j,k,Switch_No) == I) || (S_Out(j,k,Switch_No) == J_A(J))
                            m=m+1;
                        end
                    end
                end
                if m == 2     % m==2 means that the output has I and J_A(J). ˵���������м���I��Ҳ��J_A(J),��ѭ���ֵ��ĸ�����
                    %{
                    %[start monitoring last_ij, if these two flows keep existing, then it is "true", pass. ]
                    %[if they do not exist at the same time, then the last one is the "last", and test one more time for port]
                    whatchdog = 1;      %start monitoring
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
                disp(['last_',num2str(I),num2str(J_A(J)),' = S',num2str(last_ij),num2str(last_ij_port),', and the switch is S',num2str(last_ij)]);
                %store first_ij with structure Cmpt.first
                %query first_ij
                %without considering I = J_A(J)
                Cmpt.last(I,J_A(J),1) = last_ij;
                Cmpt.last(I,J_A(J),2) = last_ij_port;
                break
            end
        end
        judge_last = 0;
    end
    clear judge_last J i Switch_No j m k last_ij last_ij_port hang lie;
end
disp('********************************************');
disp(' ');


%% ***************************[Useless]*********************************
%{
%Cmpt(1,:) = cmpt;
%Cmpt(2,:) = cmpt_out;
%Cmpt
%�ֶ������������
%CMPT = input('\n��ο����˽ṹͼ�����ݴ���I��J�����Ľ�������˳�򣬶Ա�Cmpt�ĵ�һ��Ԫ����ţ��Ծ������ʽ����һ��Ϊ��λ�����Ⱥ�˳������¾���\n'); %Cmpt��һ�е������Ǽ��ʹ���S��
%CMPT = [1 2;1 1];
%CMPT1 = size(CMPT);  %ȡCMPT�Ĺ�������̶�Ϊ2������������ھ���������˿���

disp(' ');
disp('****************************');
disp(['first_ij = S',num2str(CMPT(1,1)),num2str(CMPT(2,1))]);       %CMPT(1,1)��������S����CMPT(2,1)��������˿ڼ�
disp(['last_ij = S',num2str(CMPT(1,CMPT1(2))),num2str(CMPT(2,CMPT1(2)))]);
disp('���¼��first_ij��last_ij');
disp('****************************');
disp(' ');
%}


%% ************************[Useless or not???]***********************
%{
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


%% ********************* Calculate [S_first_ij_min_j] and [S_first_ij_min_i]  still need the calculation of A **************************
%{
%show the value of first_ij firstly, like S11
%disp(['first_ij = S',num2str(CMPT(1,1)),num2str(CMPT(2,1))]);       %CMPT(1,1) means S_X, CMPT(2,1) means output port
%then find out the subscript switch of first_ij( = No_first_ij), eg. the 2nd one, we need to find the location of CMPT(1,1)
%n = 1;
%while ~isequal(S_In(:,:,CMPT(1,1)),S_InSeq(:,:,n))      %compare S[CMPT(1,1)] with the switch 
                                                         %in the topology diagram, to find the diffence and same.
%        n = n + 1;      %when find the first_ij, and the location of it(n), jump out the loop.
%end
%No_first_ij = n;
%}

%S.max = zeros(Max_Output,Tao,Num);   %save all the S pararmeters, page---I/J, number = ��, row---the name of switch
%S.min = zeros(Max_Output,Tao,Num);   %colomn---the number of output port.

%for Class B, find out the location that J go in(the location will be No_J), like the 1st.[find J first, then other elements]
%***********************************************************
%disp('This is the calculation towards to Calss B: ');
disp('*****************************');
disp('Calculating the S.min ...');
disp('*****************************');
disp(' ');
count_j = 1;
for J=1:length(J_RealB)
    for i=1:Num                             %page
        for j=1:Max_Input                   %row
            n = 0;
            N = 0;
            for k=1:Max_InputNum            %colomn
                if S_InSeq(j,k,i) == J_RealB(J)      %find J_B one by one
                    n=n+1;
                end
                if S_InSeq(j,k,i) == I      %find I one by one
                    N=N+1;
                end
            end
            %********************
            if n == 1                       %this means there is J at this row(this input port)
                m = 0;
                for l = 1:Max_InputNum
                    if (S_InSeq(j,l,i) ~= J_RealB(J)) && (S_InSeq(j,l,i) ~= 0)     %this means there is non-J at this row(this input port)
                        m=m+1;
                    end
                end
                if m == 0                   %do not exist non-J element(only J)
                    No_J_B(count_j) = i;      %No_J means that the location of switch in S_InSeq... Chinese: ��������ʾS_InSeq�еĵڼ���������
                    count_j = count_j + 1;
                end
            end
            %********************
            if N == 1                       %there is I in this row(this input port)
                M = 0;
                for L=1:Max_InputNum
                    if (S_InSeq(j,L,i) ~= I) && (S_InSeq(j,L,i) ~= 0)     %this means there is non-I at this row(this input port)
                        M=M+1;
                    end

                end
                if M == 0                %do not exist non-I element(only I)
                    No_I = i;
                end
            end
            %********************
        end
    end
end

%claculate S_first_ij_min_j
%S_first_ij_min_j = C��I��*��No_first_ij - No_J_B + 1��

%disp('**************[Calculate S_first_ij_min_j]*******************');
%disp('********************************************************');

for count_j = 1:length(J_RealB)
    %Name_first_ij = input(['Please enter the first number in order, first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'�Ľ�������źͶ˿ڣ�������������ʽ���룺\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealB(count_j),1)),S_InSeq(:,:,n))       %compare S[CMPT(1,1)] with the switch 
                                                                                    %in the topology diagram, to find the diffence and same.
        n = n + 1;      %when find the first_ij, and the location of it(n), jump out the loop.
    end
    %save the location, name and port of first_ij
    %see the parameter of S.min
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
    %disp('--------------------------');
    %disp(' ');
    %S.min(Name_first_ij(1),Name_first_ij(2),J_RealB(count_j)) =...
    %    C(J_RealB(count_j)) * (No_first_ij - No_J(count_j) + 1);
    S.min( Cmpt.first(I,J_RealB(count_j),1) , Cmpt.first(I,J_RealB(count_j),2) , J_RealB(count_j)) =...
        C(J_RealB(count_j)) * (No_first_B(I,J_RealB(count_j)) - No_J_B(count_j) + 1);
end
%disp('*********************************');
%disp(' ');
%disp(' ');

%����S_first_ij_min_i
%S_first_ij_min_i = C(I) * (No_first_ij - No_I + 1)
%disp('**************[Calculate the S_first_ij_min_i]*******************');
%disp('********************************************************');
%disp(['This is the calculation of I, and I = ',num2str(I)]);

%{
Do not find out the function.

if Pr(I) == 3
    disp('I belongs to Class A');
else if Pr(I) == 2
        disp('I belongs to Class B');
    else
        disp('I belongs to Class C');
    end
end
%}
%disp(' ');
for count_j = 1:length(J_RealB)
    %Name_first_ij = input(['�밴˳�����֮ǰ��¼��first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'�Ľ�������źͶ˿ڣ�������������ʽ���룺\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealB(count_j),1)),S_InSeq(:,:,n))      
        n = n + 1;      %when find the first_ij, and the location of it(n), jump out the loop.
    end
    No_first_B(I,J_RealB(count_j)) = n;
    %No_first_ij = input(['first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'���ڵ�S',...
    %    num2str(Cmpt.first(I,J_RealB(count_j),1)),'�ǵڼ�����������\nNo_first_ij = ']);
    %disp('�����¼����');
    %disp('--------------------------');
    %eval(['S_first_',num2str(I),num2str(J_RealB(count_j)),'_min_',num2str(I),' = ',...
    %    num2str(C(J_RealB(count_j)) * (No_first(I,J_RealB(count_j)) - No_I + 1))]);
    %disp('--------------------------');
    %disp(' ');
    %S.min(Name_first_ij(1),Name_first_ij(2),I) =...
    %    C(J_RealB(count_j)) * (No_first_ij - No_I + 1);
    %S.min( Cmpt.first(I,J_RealB(count_j),1) , Cmpt.first(I,J_RealB(count_j),2) , I) =...
    %    C(J_RealB(count_j)) * (No_first_ij - No_I + 1);     %S.min����S_mnΪ��λ�Ĵ洢�ռ䣬��������first/lastΪ��λ�Ĵ洢�ռ�
    S.min( Cmpt.first(I,J_RealB(count_j),1) , Cmpt.first(I,J_RealB(count_j),2) , I) =...
        C(J_RealB(count_j)) * (No_first_B(I,J_RealB(count_j)) - No_I + 1);     %S.min is a store space with the unit of S_mn,
                                                                               %rather than the unit of first/last
                                                                               %Chinese: S.min����S_mnΪ��λ�Ĵ洢�ռ䣬��������first/lastΪ��λ�Ĵ洢�ռ�
end
%disp('*********************************');
%disp(' ');
%disp(' ');

%***********************************************************
%disp('This is the calculation towards to Class A: ');
count_j = 1;
for J=1:length(J_RealA)
    for i=1:Num                             %page
        for j=1:Max_Input                   %row
            n = 0;
            N = 0;
            for k=1:Max_InputNum            %colomn
                if S_InSeq(j,k,i) == J_RealA(J)      %find J_A one by one
                    n=n+1;
                end
                if S_InSeq(j,k,i) == I      %find I one by one
                    N=N+1;
                end
            end
            %********************
            if n == 1                       %this means there is J at this row(this input port)
                m = 0;
                for l = 1:Max_InputNum
                    if (S_InSeq(j,l,i) ~= J_RealA(J)) && (S_InSeq(j,l,i) ~= 0)     %this means there is non-J at this row(this input port)
                        m=m+1;
                    end
                end
                if m == 0                   %do not exist non-J element(only J)
                    No_J_A(count_j) = i;      %No_J means that the location of switch in S_InSeq... Chinese: No_J��������ʾS_InSeq�еĵڼ���������
                    count_j = count_j + 1;
                end
            end
            %********************
            if N == 1                       %there is I in this row(this input port)
                M = 0;
                for L=1:Max_InputNum
                    if (S_InSeq(j,L,i) ~= I) && (S_InSeq(j,L,i) ~= 0)     %this means there is non-I at this row(this input port)
                        M=M+1;
                    end

                end
                if M == 0                %do not exist non-I element(only I)
                    No_I = i;
                end
            end
            %********************
        end
    end
end

%Calculate S_first_ij_min_j
%S_first_ij_min_j = C��I��*��No_first_ij - No_J_A + 1��
%disp('**************[Calculate the S_first_ij_min_j]*******************');
%disp('********************************************************');
for count_j = 1:length(J_RealA)
    %Name_first_ij = input(['�밴˳�����֮ǰ��¼��first_',num2str(I),...
    %    num2str(J_RealB(count_j)),'�Ľ�������źͶ˿ڣ�������������ʽ���룺\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealA(count_j),1)),S_InSeq(:,:,n))      %Compare S_In with S_InSeq
        n = n + 1;      %when find the first_ij, and the location of it(n), jump out the loop.
    end
    %save the location, name and port of first_ij
    %see the parameter of S.min
    No_first_A(I,J_RealA(count_j)) = n;     %location
    Name_first_A(I,J_RealA(count_j)) = Cmpt.first(I,J_RealA(count_j),1);    %name
    Port_first_A(I,J_RealA(count_j)) = Cmpt.first(I,J_RealA(count_j),2);    %port
    
    %disp('--------------------------');
    %disp(' ');
    S.min( Cmpt.first(I,J_RealA(count_j),1) , Cmpt.first(I,J_RealA(count_j),2) , J_RealA(count_j)) =...
        C(J_RealA(count_j)) * (No_first_A(I,J_RealA(count_j)) - No_J_A(count_j) + 1);
end
%disp('*********************************');
%disp(' ');
%disp(' ');

%Calculate S_first_ij_min_i
%S_first_ij_min_i = C(I) * (No_first_ij - No_I + 1)
%disp('**************[Calculate S_first_ij_min_i]*******************');
%disp('********************************************************');
%disp(['This is the calculation of I, and I = ',num2str(I)]);

%{
if Pr(I) == 3
    disp('I belongs to Class A');
else if Pr(I) == 2
        disp('I belongs to Class B');
    else
        disp('I belongs to Class C');
    end
end
%}
%disp(' ');
for count_j = 1:length(J_RealA)
    %Name_first_ij = input(['�밴˳�����֮ǰ��¼��first_',num2str(I),...
    %    num2str(J_RealA(count_j)),'�Ľ�������źͶ˿ڣ�������������ʽ���룺\n']);
    n = 1;
    while ~isequal(S_In(:,:,Cmpt.first(I,J_RealA(count_j),1)),S_InSeq(:,:,n))     
        n = n + 1;      %when find the first_ij, and the location of it(n), jump out the loop.
    end
    No_first_A(I,J_RealA(count_j)) = n;
    %disp('--------------------------');
    %disp(' ');
    S.min( Cmpt.first(I,J_RealA(count_j),1) , Cmpt.first(I,J_RealA(count_j),2) , I) =...
        C(J_RealA(count_j)) * (No_first_A(I,J_RealA(count_j)) - No_I + 1);     %S.min is a store space with the unit of S_mn,
                                                                               %rather than the unit of first/last
                                                                               %Chinese: S.min����S_mnΪ��λ�Ĵ洢�ռ䣬��������first/lastΪ��λ�Ĵ洢�ռ�
end
%disp('*********************************');
%disp(' ');
%disp(' ');
clear i j k l L m M n N count_j J;


%% ***********************������S_first_ij_max_j��S_first_ij_max_i��*******************
disp('*****************************');
disp('Calculating the S.max ...');
disp('*****************************');
disp(' ');
 %{

%�������ƣ�Ӧ������Pi����
CreditA = 0;
CreditB = 0;
[hang lie] = size(Pi);  %hangû��ʲô���壬lie����һ���������ϵ�����˿ڣ��Ⱦ���lie���پ���hang��
%for i = 1:lie
    i = 1;
    S_research_temp = S_Out(Pi(2,i),:,Pi(1,i));  %����S���������������˿ڵ�������ݣ����Ǹ�һά����  [2 0 3 6 0]
    S_research_local = find(S_research_temp ~= 0);  %����ĳ������˿���˵S_research_local = [1 3 4]
    for localS = 1:length(S_research_local)         %
        %S_research = [2 3 6]  ��0Ԫ�ص��������
        
        S_research(localS) = S_Out(Pi(2,i),S_research_local(localS),Pi(1,i));   %localS = 2ʱ��S_research_local = 3
    end

        %clear S_research_temp S_research_local localS
        %*****************************
        %��ʼ��Class
        ClassA = [];
        ClassB = [];
        ClassC = [];
        for k = 1:length(S_research)
            Pr_research(k) = Pr(S_research(k));
        end
        %for k = 1:length(S_research)
        %    Pr_research(k) = Pr(S_research(k)); %��ø�����˿��������ݵ����ȼ�
            ClassA_temp = find(Pr_research == 3);   %��ø�����˿������ȼ�=3�ĽǱ�λ��
            for localA = 1:length(ClassA_temp)
            ClassA(localA) = S_Out(Pi(2,i),ClassA_temp(localA),Pi(1,i));    %���ClassA������
            end
        %end
        %for k = 1:length(S_research)
        %    Pr_research(k) = Pr(S_research(k)); %��ø�����˿��������ݵ����ȼ�
            ClassB_temp = find(Pr_research == 2);   %��ø�����˿������ȼ�=2�ĽǱ�λ��
            for localB = 1:length(ClassB_temp)
                ClassB(localB) = S_Out(Pi(2,i),ClassB_temp(localB),Pi(1,i));
            end
        %end
        %for k = 1:length(S_research)
        %    Pr_research(k) = Pr(S_research(k)); %��ø�����˿��������ݵ����ȼ�
            ClassC_temp = find(Pr_research == 1);   %��ø�����˿������ȼ�=1�ĽǱ�λ��
            for localC = 1:length(ClassC_temp)
                ClassC(localC) = S_Out(Pi(2,i),ClassC_temp(localC),Pi(1,i));
            end
        %end
        %clear ClassA_temp ClassB_temp ClassC_temp localA localB localC k
        %*****************************
        %��S_research�еĵ�j��Ԫ�طŵ�������Class�����е����һ��
    for j = 1:length(S_research)
        delaytime = 0;
        if Pr(S_research(j)) == 3
            a = find(ClassA == S_research(j));  %S_research�еĵ�j��Ԫ����ClassA�е�λ��
            temp = ClassA(a);
            ClassA(a) = []; %ɾ�����Ԫ��
            ClassA = [ClassA temp]; %�������j��Ԫ���õ����
        else if Pr(S_research(j)) == 2
                b = find(ClassB == S_research(j));  %S_research�еĵ�j��Ԫ����ClassA�е�λ��
                temp = ClassB(b);
                ClassB(b) = []; %ɾ�����Ԫ��
                ClassB = [ClassB temp]; %�������j��Ԫ���õ����
             else if Pr(S_research(j)) == 1
                c = find(ClassC == S_research(j));  %S_research�еĵ�j��Ԫ����ClassA�е�λ��
                temp = ClassC(c);
                ClassC(c) = []; %ɾ�����Ԫ��
                ClassC = [ClassC temp]; %�������j��Ԫ���õ����
                end
            end
        end
        %clear a b c temp
        %*****************************
        %ȡ�ղ�����Class������Ե�C��Pr
        for tiqu = 1:length(ClassA) %tiqu��ζ����ȡ
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
        %��ʽ��ʼ����
        %Ŀǰֻ����ֻ��B�����
        if ~isempty(ClassA)     %���ClassA���鲻�ǿռ�
            disp('��ʱ�������������');   
        else if ~isempty(ClassB)        %���ClassA�����ǿռ� �� ClassB���鲻�ǿռ�
                cc = 1;     %��ClassC�м���
                for trans = 1:(length(ClassB)-1)
                    if CreditB>=0
                        %{
                        if ~isempty(ClassC) && trans == length(ClassB)  %�ӱ��о����������еĵ������
                            delaytime = delaytime + C(ClassC(cc));
                            cc = cc + 1;
                            CreditB = CreditB + C(ClassC(cc)) * alfa_pos_B;
                        end
                        %}
                        delaytime = delaytime + C(ClassB(trans));   %����ClassB��֡
                        CreditB = CreditB + C(ClassB(trans)) * alfa_neg_B; %������ɺ��CreditB
                    end
                    if CreditB<0
                        Replenish_B = (0 - CreditB) / (alfa_pos_B); %�Ӹ�ֵ�ָ���0����ʱ
                        if ~isempty(ClassC) %&& (Replenish_B >= C(ClassC(cc))) ���C�ǿ�

                            delaytime = delaytime + Replenish_B;    %CreditB�ڻָ�ͬʱ����������C
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
        
        %�ִ����⣺�����1��2��3��һ��������������������������
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
%S.max(3,1,3) = 40;%��һ���ԣ���I=1��������һ�У�
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

%}


%% ***********************[M_h_i]*************************
disp('*****************************');
disp('Calculating the M_h_i ...');
disp('*****************************');
disp(' ');
%�����ˣ�����MӦ�ð���first_ij���㡿
%{
%����MӦ����Ҫ����Pi����  <----  �Ǵ���뷨
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
for count_m = 1:length(J_RealA)
    %m_h_i�ǲ�����ECU�����翪ʼʱ��
    %�洢������M_first_ij_i����ʽ������������M_S_mn_i����ʽ
    M_h_i_A(Cmpt.first(I,J_RealA(count_m),1) , Cmpt.first(I,J_RealA(count_m),2) , I) =...
        (S_SwitchSeq(Cmpt.first(I,J_RealA(count_m),1)) - Flow_Source(I) + 1) * C(I);
end
M_ECU_i = 0;

clear count_m 


%% ******************[All initiate set down]******************
disp('**********************************');
disp('All initiate parameters set down');
disp('**********************************');
disp(' ');


%% **************************��A_ij��**************************
%����Class B��A_ij
for i = 1:length(J_RealB)
    % A_ij = S_first_ij_max_j - M_first_ij_i
    % Name_first(I,J_RealB(i))�����˽�������ţ����ƣ���Ҳ����S_mn�е�m
    % Port_first(I,J_RealB(i))����������˿ڣ�Ҳ����S_mn�е�n
    A(I,J_RealB(i)) = S.max(Name_first_B(I,J_RealB(i)) , Port_first_B(I,J_RealB(i)) , J_RealB(i)) -...
        M_h_i_B(Name_first_B(I,J_RealB(i)) , Port_first_B(I,J_RealB(i)) , I);
end
%����Class A��A_ij
for i = 1:length(J_RealA)
    A(I,J_RealA(i)) = S.max(Name_first_A(I,J_RealA(i)) , Port_first_A(I,J_RealA(i)) , J_RealA(i)) - ...
        M_h_i_A(Name_first_A(I,J_RealA(i)) , Port_first_A(I,J_RealA(i)) , I);
end


%% ***********************[D_ClassX]************************
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


%% *****************************[Rep_i_X]*******************************
if ismember(I,J_B)
    Rep_B(I) = Sum_D_ClassB * (alfa_neg_B / alfa_pos_B);
else if ismember(I,J_A)
        Rep_A(I) = Sum_D_ClassA * (alfa_neg_A / alfa_pos_A);
    end
end


%% ************************[���е��������Ĳ���]********************************
%���¶���hp���ԣ�ֻ��X == B��ʱ��ų�����������ʵ���ڼ���D_hp_B
if ismember(I,J_B)
    %Ҫ��Ū�����е�last_ij��j����hpX����Class A����Ϊ�˴�������I��Class C��
    Temp = 0;   %��ֵ���ΪD_hp_B
    for i = 1:length(J_A)
        disp(['last_',num2str(I),num2str(J_A(i)) , '��S' , num2str(Cmpt.last(I,J_A(i),1)) ,...
            num2str(Cmpt.last(I,J_A(i),2))]);
        eval(['syms W_S',num2str(Cmpt.last(I,J_A(i),1)),num2str(Cmpt.last(I,J_A(i),2)),...
            '_',num2str(I)]);
        %��֪��Ҫ��ô����
        %eval(['temp = positive(1 + floor( (W_S',num2str(Cmpt.last(I,J_A(i),1)),num2str(Cmpt.last(I,J_A(i),2)),...
        %    '_',num2str(I),') / T(J_A(i)) )) * C(J_A(i))']);
    end
end


%% ***********************���Ӻ�����********************************
%�ⲿ�������ϴ�����ʹ�ù����Ӻ���
%*********************************************************
% I - TRANST


%*********************************************************





%% �� D_lpX_i


%% �� max(Cj)


%% �� (num(Pi)-1)*sl


%% �� delta_h_i_t


%% �� max(alfa_neg_X / alfa_pos_X)

