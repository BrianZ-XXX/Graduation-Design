clear all;
clc;
% ��
% �ڲ�������m�ļ�֮ǰ�����ò��������д˳��򡣣���λ��us��
% �ٴ��ļ��е�ȫ�ֱ���I���޸Ĳ����������(ֱ��������global I)

%% չʾ�������˽ṹ
%map=imread('topology_eg.jpg');
%imshow(map);

%% ����ģ�Ͷ�����֪������ 
disp('ע�⣡');
disp('�������������Ӧ��Ϊ�����������ʽ��');
global C;
       %C=[40 40 40 40 40 120];                 %����ʱ�䣨����֡����
       C=input('\n������ÿһ���������ĵ�λ֡����ʱ�䣺\n');
global T;
       %T=[2000 2000 2000 160 2000 2000];       %������
       T=input('\n������ÿһ���������ĵ�λ֡������ʱ�䣺\n');
global Pr;
       %Pr=[2 2 2 3 1 1];                       %���ȼ�
       Pr=input('\n��Ӧ�1����n��˳������ÿһ�������������ȼ���\n');
% Ӧ�þ�������������Ŀ�� 
%*************************************************************
%S������
%{
%��SŪ����ά����
global S;
       S=zeros(3,3,3);  %�д�����༸������˿ڣ��д���ĳ������˿���������������ҳ������������
       S(1,:,1)=[1 0 0];
       S(2,:,1)=[2 3 6];
       S(1,:,2)=[1 2 3];
       S(2,:,2)=[4 0 0];
       S(3,:,2)=[5 0 0];
       S(1,:,3)=[2 0 0];
       S(2,:,3)=[3 0 0];
       S(3,:,3)=[6 0 0];
       %S(:,:,1)=input('\n������S1��\n��ȷ����������Ǿ����������ʽ��\n�磺[a b;c d]  ÿһ�д���һ������˿ڡ�\n');
       %S(:,:,2)=input('\n������S2��\n��ȷ����������Ǿ����������ʽ��\n�磺[a b;c d]  ÿһ�д���һ������˿ڡ�\n');
       %S(:,:,3)=input('\n������S3��\n��ȷ����������Ǿ����������ʽ��\n�磺[a b;c d]  ÿһ�д���һ������˿ڡ�\n');
       S(:,:,1)
       S(:,:,2)
       S(:,:,3)

global S1;                                      %������������
      %S1=[1 0 0;2 3 6];
       S1=input('\n������S1��\n��ȷ����������Ǿ����������ʽ��\n�磺[a b;c d]  ÿһ�д���һ������˿ڡ�\n');
global S2;
       %S2=[1 2 3;4 0 0;5 0 0];
       S2=input('\n������S2��\n��ȷ����������Ǿ����������ʽ��\n�磺[a b;c d]  ÿһ�д���һ������˿ڡ�\n');
global S3;
       %S3=[2;3;6];
       S3=input('\n������S3��\n��ȷ����������Ǿ����������ʽ��\n�磺[a b;c d]  ÿһ�д���һ������˿ڡ�\n');
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
%global F;              %�ƺ����� 
%       F=[4000 4000 4000 4000 4000 12000];
global S31;
global S11;
global S21;

global sl;          %switch latency  (����������ʹ��)
       sl=0;
       
global I;
       I=2;         %��ӳٶ���ģ����ȼ����Ǳꣿ(������ļ�Line 5)

global SeriaTerm;
global S_S21_min4;
global A_1_4;
global D_classB_1_j_t;
global Rep_1_B;
global delta_S21_1

t=0;
C
%% I - �������ģ��
S_OUT1=Sequence(S1);        %S_OUT1����['S1���������';pr].  ������һ�������������������һ��Ԫ��Ϊ��S1��������У�
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

%% II - ��h�˿ڦ�_i��������ܿ�ʼʱ��
% ��һ�����ǹ��ڼ���ģ�͵ģ�ֻ��Ҫ���������� 

%M_h_i����·��P_i��h�˿ڴ��Ħ�_i֡��busy period��������ܿ�ʼʱ�䣬��ζ�š�����ʱ�䡱��
disp('--------------------------------------');
disp('ע�⣡   M_h_i������·��P_i��h�˿ڴ��Ħ�_i֡��');
disp('        busy period��������ܿ�ʼʱ�䣬��ζ�š�����ʱ�䡱��');
M_ECU1_1=input('\n ������M_ECU1_1��\n �����ڴ����˽ṹ���ֵΪ0����\n');
M_S11_1=input('\n Please Enter M_S11_1, \nfor this topology the value is C(1): \n');
M_S21_1=input('\n Please Enter M_S21_1, \nfor this topology the value is 2*C(1): \n');

%% III - ��first_i����һ���������h�����/��С�ӳ�
% ��һ�����ǹ��ڼ���ģ�͵ģ�ֻ��Ҫ��������

%S_h_maxi/S_h_mini����������_i�е�һ��֡��frist_iһֱ��������˿�h�����/��С�ӳ١�
%first_i����������_i������·��P_i�������ĵ�һ������˿ڡ�
disp(' ');
disp('-------------------------------------');
disp('Attention: S_h_maxi/S_h_mini: Maximum/Minimum delay of a frame of flow');
disp('      ��_i from first_i till its arrival at h');
disp('first_i��First output port visited by flow ��_i along its path P_i');
S_S11_min1=input('\nPlease Enter S_S11_min1, \nfor this topology the value is C(1): \n');
S_S11_max1=input('\nPlease Enter S_S11_max1, \nfor this topology the value is C(1): \n');
S_S11_min2=input('\nPlease Enter S_S11_min2, \nfor this topology the value is 2*C(2): \n');

%% IV - Decouple the Model And Revise Any Parameter If Any Operation Is Needed
%S_S11_max2=CBS(S11,pr1)+C(S11(length(S11)));
S_S11_max2=input('\nPlease Enter S_S11_max2, the format is like: CBS(S?1,pr?)��\nPlease Enter: ') + C(S11(length(S11)));         %you can input the sequence like CBS(S11,pr1) directly...
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

%% VI - delay introduced by flow ��4(hpB)
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

%% XII - the delay of flow ��1
R1=W_S21_1_t + C(1) - t
pessimism=R1 - (W_S21_1_t + 0)      %0 is the real I_S11_Binit

%% Test Part
%global Pr31;       %Sequence, Line 14

%% End