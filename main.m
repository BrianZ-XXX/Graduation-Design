clear all;
clc;
% τ
% 在测试其他m文件之前先设置参数并运行此程序。（单位：us）
% 再此文件中的全局变量I中修改并决定最坏对象。(直接搜索：global I)

%% 展示网络拓扑结构
%map=imread('topology_eg.jpg');
%imshow(map);

%% 根据模型定义已知参数。 
disp('注意！');
disp('下面输入的数据应该为矩阵或向量形式！');
global C;
       %C=[40 40 40 40 40 120];                 %传输时间（代表帧长）
       C=input('\n请输入每一组数据流的单位帧传输时间：\n');
global T;
       %T=[2000 2000 2000 160 2000 2000];       %传输间隔
       T=input('\n请输入每一组数据流的单位帧传输间隔时间：\n');
global Pr;
       %Pr=[2 2 2 3 1 1];                       %优先级
       Pr=input('\n请从τ1到τn的顺序输入每一组数据流的优先级：\n');
% 应该决定交换机的数目。 
%*************************************************************
%S的输入
%{
%把S弄成三维矩阵
global S;
       S=zeros(3,3,3);  %行代表最多几个输入端口，列代表某个输入端口最多的输入数量，页代表交换机数量
       S(1,:,1)=[1 0 0];
       S(2,:,1)=[2 3 6];
       S(1,:,2)=[1 2 3];
       S(2,:,2)=[4 0 0];
       S(3,:,2)=[5 0 0];
       S(1,:,3)=[2 0 0];
       S(2,:,3)=[3 0 0];
       S(3,:,3)=[6 0 0];
       %S(:,:,1)=input('\n请输入S1，\n并确定你的输入是矩阵或向量形式，\n如：[a b;c d]  每一行代表一个输入端口。\n');
       %S(:,:,2)=input('\n请输入S2，\n并确定你的输入是矩阵或向量形式，\n如：[a b;c d]  每一行代表一个输入端口。\n');
       %S(:,:,3)=input('\n请输入S3，\n并确定你的输入是矩阵或向量形式，\n如：[a b;c d]  每一行代表一个输入端口。\n');
       S(:,:,1)
       S(:,:,2)
       S(:,:,3)

global S1;                                      %设置输入序列
      %S1=[1 0 0;2 3 6];
       S1=input('\n请输入S1，\n并确定你的输入是矩阵或向量形式，\n如：[a b;c d]  每一行代表一个输入端口。\n');
global S2;
       %S2=[1 2 3;4 0 0;5 0 0];
       S2=input('\n请输入S2，\n并确定你的输入是矩阵或向量形式，\n如：[a b;c d]  每一行代表一个输入端口。\n');
global S3;
       %S3=[2;3;6];
       S3=input('\n请输入S3，\n并确定你的输入是矩阵或向量形式，\n如：[a b;c d]  每一行代表一个输入端口。\n');
%}          

%% *************************************************************
global aA_pos;
       aA_pos=50;
global aA_neg;
       aA_neg=-50;
global aB_pos;
       aB_pos=25;
global aB_neg;
       aB_neg=-75;
%global F;              %似乎无用 
%       F=[4000 4000 4000 4000 4000 12000];
global S31;
global S11;
global S21;

global sl;          %switch latency  (其他函数中使用)
       sl=0;
       
global I;
       I=2;         %最坏延迟对象的？优先级？角标？(请见本文件Line 5)

global SeriaTerm;
global S_S21_min4;
global A_1_4;
global D_classB_1_j_t;
global Rep_1_B;
global delta_S21_1

t=0;
C
%% I - 输出序列模型
S_OUT1=Sequence(S1);        %S_OUT1代表['S1的输出序列';pr].  （这是一个两行向量，并且其第一行元素为：S1的输出序列）
S_OUT2=Sequence(S2);
S_OUT3=Sequence(S3);
Num1=size(S_OUT1);
for i=1:Num1(2)
    S11(i)=S_OUT1(2*i-1);
    pr1(i)=Pr(S11(i));
end
Num2=size(S_OUT2);
for i=1:Num1(2)
    S21(i)=S_OUT2(2*i-1);
    pr2(i)=Pr(S21(i));
end
Num3=size(S_OUT3);
for i=1:Num3(2)
    S31(i)=S_OUT3(2*i-1);
    pr3(i)=Pr(S31(i));
end

%S11=S_OUT1(1);
%S21=S_OUT2(1);
%S31=S_OUT3(1);

%% II - 在h端口τ_i的最早可能开始时间
% 这一部分是关于计算模型的，只需要调整参数。 

%M_h_i：在路径P_i中h端口处的τ_i帧的busy period的最早可能开始时间，意味着“最少时间”。
disp('--------------------------------------');
disp('注意！   M_h_i代表：在路径P_i中h端口处的τ_i帧的');
disp('        busy period的最早可能开始时间，意味着“最少时间”。');
M_ECU1_1=input('\n 请输入M_ECU1_1，\n （对于此拓扑结构这个值为0）：\n');
M_S11_1=input('\n Please Enter M_S11_1, \nfor this topology the value is C(1): \n');
M_S21_1=input('\n Please Enter M_S21_1, \nfor this topology the value is 2*C(1): \n');

%% III - 从first_i（第一个输出）到h的最大/最小延迟
% 这一部分是关于计算模型的，只需要调整参数

%S_h_maxi/S_h_mini：数据流τ_i中的一个帧从frist_i一直到它到达端口h的最大/最小延迟。
%first_i：数据流τ_i沿着其路径P_i所经历的第一个输出端口。
disp(' ');
disp('-------------------------------------');
disp('Attention: S_h_maxi/S_h_mini: Maximum/Minimum delay of a frame of flow');
disp('      τ_i from first_i till its arrival at h');
disp('first_i：First output port visited by flow τ_i along its path P_i');
S_S11_min1=input('\nPlease Enter S_S11_min1, \nfor this topology the value is C(1): \n');
S_S11_max1=input('\nPlease Enter S_S11_max1, \nfor this topology the value is C(1): \n');
S_S11_min2=input('\nPlease Enter S_S11_min2, \nfor this topology the value is 2*C(2): \n');

%% IV - Decouple the Model And Revise Any Parameter If Any Operation Is Needed
%S_S11_max2=CBS(S11,pr1)+C(S11(length(S11)));
S_S11_max2=input('\nPlease Enter S_S11_max2, the format is like: CBS(S?1,pr?)，\nPlease Enter: ') + C(S11(length(S11)));         %you can input the sequence like CBS(S11,pr1) directly...
S_S21_max4=input('\nPlease Enter S_S21_max4, \nfor this topology the value is C(4): \n');
S_S21_min4=input('\nPlease Enter S_S21_min4, \nfor this topology the value is C(4): \n');

%% V - the workload introduced by all classB flows
A_1_1 = 0;
A_1_2 = S_S11_max2 - M_S11_1;
A_1_3 = A_1_2;
D_classB_1_1_t = positive(1+floor((t + 0 - 0 + A_1_1) / T(1))) * C(1);
D_classB_1_2_t = positive(1+floor((t + S_S11_max1 - S_S11_min2 + A_1_2) / T(2))) * C(2);
D_classB_1_3_t = positive(1+floor((t + S_S11_max1 - S_S11_min2 + A_1_3) / T(3))) * C(3);
D_classB_1_j_t = D_classB_1_1_t + D_classB_1_2_t + D_classB_1_3_t;

%% VI - delay introduced by flow τ4(hpB)
A_1_4=S_S21_max4 - M_S21_1;
%D_hpB_1_t = positive(1+floor((W_S21_1_t-S_S21_min4+A_1_4)/T(4))) * C(4);
%D_Amax_1_t=D_hpB_1_t;
%D_Amin_1_t=min(D_hpB_1_t , ((W_S21_1_t-S_S21_min4+A_1_4)*(aA_pos/(aA_pos+abs(aA_neg)))));

%% VII - the maximum replenshing time of classB flows
Rep_1_B=D_classB_1_j_t * (abs(aB_neg)/aB_pos);

%% VIII - the maximum workload intruduced by all classB flows
%D_Bmax_1_t=D_classB_1_j_t+positive(Rep_1_B-D_Amin_1_t);

%% IX - the maximum value of l_S11_0 and the minumum value of l_S11_1
k=1;
for i=1:length(C)
    if (Pr(i)==2)||(Pr(i)==3)
        cc=C(k);
        k=k+1;
    end
end
C_IPs11_0_Amin=min(cc);

k=1;
for i=1:length(C)
    if Pr(i)==2
        cc=C(k);
        k=k+1;
    end
end
C_IPs11_1_Bmax=max(cc);

max_l_S11_0=(positive(1 + floor(t/T(4))) * C(4) - C_IPs11_0_Amin) * (1+abs(aB_neg)/aB_pos);
min_l_S11_1=positive(positive(1+floor((t + S_S11_max1 - S_S11_min2 + A_1_2) / T(2))) * C(2) + ...
    positive(1 + floor((t + S_S11_max1 - S_S11_min2 + A_1_3) / T(3))) * C(3) - C_IPs11_1_Bmax);

%% X - serialization term
delta_S11_1=0;
delta_S21_1=C(5);
k=1;
for i=1:length(C)
    if Pr(i)==2
        cc=C(k);
        k=k+1;
    end
end
C_S11_Bmax=max(cc);

I_S11_Binit=C_S11_Bmax*(abs(aA_neg)/aA_pos);
Delta_S11_1_t=positive(min_l_S11_1 - max_l_S11_0 - max(delta_S11_1 , I_S11_Binit));
Delta_S21_1_t=0;
SeriaTerm=positive(Delta_S11_1_t + Delta_S21_1_t - t);

%% XI - Calculate the End-to-end Delay
W_S21_1_t=Iteration(C(2))
D_hpB_1_t = positive(1 + floor((W_S21_1_t - S_S21_min4 + A_1_4) / T(4))) * C(4); 
D_Amax_1_t=D_hpB_1_t
D_Amin_1_t=min(D_hpB_1_t , ((W_S21_1_t - S_S21_min4 + A_1_4) * (aA_pos/(aA_pos + abs(aA_neg)))));
D_Bmax_1_t=D_classB_1_j_t + positive(Rep_1_B - D_Amin_1_t)

%% XII - the delay of flow τ1
R1=W_S21_1_t + C(1) - t
pessimism=R1 - (W_S21_1_t + 0)      %0 is the real I_S11_Binit

%% Test Part
%global Pr31;       %Sequence, Line 14

%% End