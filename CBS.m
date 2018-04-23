%% Main Function: CBS
function DT=CBS(S,pr)    %���SӦ����ĳ��Switch������������У��ȸ�Ϊ�������룩��S31 S11 S21
% ����˼·����������������֡���ȼ�����A����B�������ϻ����¼���Aa��Ab��ֵ��
% ����Ŀ�ģ��õ��ӳ�ʱ�䡣 

%********************************************************
% ȫ�ֱ���
%S=[6 1 3 2];
%pr=[1 2 2 2];
global transT;
%global Pr;
    transT=TRANST(S);    %����������������[6 4 1 3 2],transT=[120 40 40 40 ]   for the unordered input sequence [6 4 1 3 2], transT=[120 40 40 40 40]
global C;
global aA_pos;
%global aA_neg;
global aB_pos;
%global aB_neg;
global CreditA;
global CreditB;
global A;
global B;
global BE;
global pprr;
       pprr=pr;

%*********************************************************
% ��ʼֵ
CreditA=0;
CreditB=0;
dt=0;
TIME=0;         % ��ʱ��

%*********************************************************
% ������������������A��B��C�����

%���в�ֹһ��C��֡ʱ��
    %���ȴ���C��֡�� 
    %ÿ����C��֡������Ϻ󣬼���CreditA�����CreditA>=0������һ�δ���A��֡��
    %A��֡ÿ������һ�Σ�����CA�����CA>=0������һ�δ���A��֡��
                             %���CA<0������CreditB�����CB>=0������һ�δ���B��֡��
                                                   %���CB<0������C��֡��
    %B��֡ÿ������һ�Σ�����CA�����CA>=0������һ�δ���A��֡��
                             %���CA<0������CreditB�����CB>=0������һ�δ���B��֡��
                                                   %���CB<0������C��֡��

%*********************************************************
% �ռ�����������и����ȼ��ĸ���
%S�Ǵ������룬���Ҷ�������������� [6 4 1 3 2]��0���Եõ���a=1, b=3, c=1
%a=0;
%b=0;
c=0;
for i=1:length(S)
    if pr(i)==1
        c=c+1;
    end
%    elseif pr(i)==2
%        b=b+1;
%    else
%        c=c+1;
%    end
end

A=0;
B=0;
BE=0;
k=1;                    %��ʱ����֪�����������������У�A B C���ȼ��ֱ�����Щλ��
m=1;
n=1;
for i=1:length(S)
    if pr(i)==3
        A(k)=i;         %���������������� [6 1 3 2], A=[]
        k=k+1;
    elseif pr(i)==2
        B(m)=i;         %���������������� [6 1 3 2], B=[2 3 4]
        m=m+1;
    else
        BE(n)=i;        %���������������� [6 1 3 2], BE=[1]
        n=n+1;
    end
end

%*********************************************************
% ���俪ʼ
aa=1;       %�����䵽��Ӧ�����ĵڼ���
bb=1;
cc=0;

%*********************************************************
% ��C��֡Ϊ�ǿ�ʱ����c~=0����ô 
if c~=0
    % transmit C;
    cc=cc+1;
    TIME=TIME+C(S(BE(cc)));
    CreditA=CreditA + aA_pos * transT(BE(cc));
    CreditB=CreditB + aB_pos * transT(BE(cc));
    
    %ѭ����ʼ 
    %�ܹ��ֻصĴ���Ϊ�����г���-1
    for i=1:length(S)-1
        Who=check(aa,bb,cc);           %�����ж�Credit��������жϴ�ʱ�ô��ݵ���˭
        if pr(Who)==3
            aa=aa+1;
        elseif pr(Who)==2
            bb=bb+1;
        else
            cc=cc+1; 
        end
        TIME=TIME+trans(Who);     %������һ���õ��Ĵ������Ĵ���ʱ�����  ---ʱ���ۼ�
        %����һ���ֻؽ���
    end
else
    for i=1:length(S)
        Who=check(aa,bb,cc);           %�����ж�Credit��������жϴ�ʱ�ô��ݵ���˭
        if pr(Who)==3
            aa=aa+1;
        elseif pr(Who)==2
            bb=bb+1;
        else
            cc=cc+1; 
        end
        TIME=TIME+trans(Who);     %������һ���õ��Ĵ������Ĵ���ʱ�����  ---ʱ���ۼ�
        %����һ���ֻؽ���
    end
end

DT=TIME;
%if length(BE)~=0
%    dt=dt+C(BE(1));
%    time=time+1;
%    CreditA=CreditA+aA_pos*transT(BE(1));
%    CreditB=CreditB+aB_pos*transT(BE(1));
%    
%    if CreditA>0
%        transmit A
%    end
end

%% �Ӻ���
%*********************************************************
% I - TRANST
function t=TRANST(S)
% ĳ��Switch������������������ʱ����

%*************************
% ȫ�ֱ���
global C;

%*************************
% ���ܲ���
for i=1:length(S)
    transt(i)=C(S(i));      %S��ÿ���������������ȫ����ʱ��
end
    
t=transt;
end

%*********************************************************
% II - check
function who=check(aa,bb,cc)
% ���������Ϊ�˼���Creditֵ������ΨһĿ�����ж���һ���ô���ABC�е�˭��

%*************************
% ȫ�ֱ���
global CreditA;
global CreditB;
global A;
global B;
global BE;

%*************************
% ���ܲ���

if A~=0
    if (CreditA>=0)
        W=A(aa);
    else
        if (CreditB>=0)
            W=B(bb);
        else
            W=BE(cc);
        end
    end
elseif B~=0
    if (CreditB>=0)
        W=B(bb);
    else
        W=BE(cc);
    end
end

who=W;              %�������who���õ������������ĽǱ���� 
end

%*********************************************************
% III - trans
function tt=trans(who)      % ����who����Ϊ�����������ĽǱ������

% transֻ��Ϊ���亯������                            "trans" works only as a transmission function. 
% ��������ѡ����ʱ�䣬�����Լ���credit�ı仯��      It not only selects the transmission time, but also calculates the change of credit. 
% Ȼ�����Ϳ��Ը���credit���ֳ�������Ķ���          Then it will distinguish the transmitted object according to Credit. 

%���в�ֹһ��C��֡ʱ��
    %���ȴ���C��֡�� 
    %ÿ����C��֡������Ϻ󣬼���CreditA�����CreditA>=0������һ�δ���A��֡��
    %A��֡ÿ������һ�Σ�����CA�����CA>=0������һ�δ���A��֡��
                             %���CA<0������CreditB�����CB>=0������һ�δ���B��֡��
                                                   %���CB<0������C��֡��
    %B��֡ÿ������һ�Σ�����CA�����CA>=0������һ�δ���A��֡��
                             %���CA<0������CreditB�����CB>=0������һ�δ���B��֡��
                                                   %���CB<0������C��֡��
                                                   
%*************************
% ȫ�ֱ���
global CreditA;
global CreditB;
global aA_pos;
global aA_neg;
global aB_pos;
global aB_neg;
global C;
global pprr;

%*************************
% ���ܲ���
delaytime=0;
if pprr(who) == 3
    delaytime = delaytime + C(who);
    CreditA = CreditA + aA_neg * C(who);
    CreditB = CreditB + aB_pos * C(who);
elseif pprr(who) == 2
    delaytime = delaytime + C(who);
    CreditA = CreditA + aA_pos * C(who);
    CreditB = CreditB + aB_neg * C(who);
elseif pprr(who) == 1
    delaytime = delaytime + C(who);
    CreditA = CreditA + aA_pos * C(who);
    CreditB = CreditB + aB_pos * C(who);
end

tt=delaytime;
end