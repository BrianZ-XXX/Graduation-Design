%% Main Function: CBS
function DT=CBS(S,pr)    %这个S应该是某个Switch的有序输出序列（先改为串行输入），S31 S11 S21
% 基本思路：如果被传输的数据帧优先级不是A或者B，则向上或向下计算Aa或Ab的值，
% 最终目的：得到延迟时间。 

%********************************************************
% 全局变量
%S=[6 1 3 2];
%pr=[1 2 2 2];
global transT;
%global Pr;
    transT=TRANST(S);    %对于无序输入序列[6 4 1 3 2],transT=[120 40 40 40 ]   for the unordered input sequence [6 4 1 3 2], transT=[120 40 40 40 40]
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
% 初始值
CreditA=0;
CreditB=0;
dt=0;
TIME=0;         % 总时间

%*********************************************************
% 对于数据流包含所有A、B、C的情况

%当有不止一个C级帧时：
    %首先传输C级帧。 
    %每次在C级帧传输完毕后，计算CreditA。如果CreditA>=0，则下一次传输A级帧；
    %A级帧每传输完一次，计算CA，如果CA>=0，则下一次传输A级帧；
                             %如果CA<0，计算CreditB：如果CB>=0，则下一次传输B级帧；
                                                   %如果CB<0，则传输C级帧。
    %B级帧每传输完一次，计算CA，如果CA>=0，则下一次传输A级帧；
                             %如果CA<0，计算CreditB：如果CB>=0，则下一次传输B级帧；
                                                   %如果CB<0，则传输C级帧。

%*********************************************************
% 收集输出数据流中各优先级的个数
%S是串行输入，并且对于无序输出序列 [6 4 1 3 2]，0可以得到：a=1, b=3, c=1
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
k=1;                    %这时可以知道在无序输入序列中，A B C优先级分别在哪些位置
m=1;
n=1;
for i=1:length(S)
    if pr(i)==3
        A(k)=i;         %对于无序输入序列 [6 1 3 2], A=[]
        k=k+1;
    elseif pr(i)==2
        B(m)=i;         %对于无序输入序列 [6 1 3 2], B=[2 3 4]
        m=m+1;
    else
        BE(n)=i;        %对于无序输入序列 [6 1 3 2], BE=[1]
        n=n+1;
    end
end

%*********************************************************
% 传输开始
aa=1;       %代表传输到相应类别里的第几个
bb=1;
cc=0;

%*********************************************************
% 当C级帧为非空时，即c~=0，那么 
if c~=0
    % transmit C;
    cc=cc+1;
    TIME=TIME+C(S(BE(cc)));
    CreditA=CreditA + aA_pos * transT(BE(cc));
    CreditB=CreditB + aB_pos * transT(BE(cc));
    
    %循环开始 
    %总共轮回的次数为，序列长度-1
    for i=1:length(S)-1
        Who=check(aa,bb,cc);           %根据判断Credit的情况来判断此时该传递的是谁
        if pr(Who)==3
            aa=aa+1;
        elseif pr(Who)==2
            bb=bb+1;
        else
            cc=cc+1; 
        end
        TIME=TIME+trans(Who);     %进行上一步得到的传输对象的传输时间计算  ---时间累加
        %至此一个轮回结束
    end
else
    for i=1:length(S)
        Who=check(aa,bb,cc);           %根据判断Credit的情况来判断此时该传递的是谁
        if pr(Who)==3
            aa=aa+1;
        elseif pr(Who)==2
            bb=bb+1;
        else
            cc=cc+1; 
        end
        TIME=TIME+trans(Who);     %进行上一步得到的传输对象的传输时间计算  ---时间累加
        %至此一个轮回结束
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

%% 子函数
%*********************************************************
% I - TRANST
function t=TRANST(S)
% 某个Switch输出的数据流，传输耗时向量

%*************************
% 全局变量
global C;

%*************************
% 功能部分
for i=1:length(S)
    transt(i)=C(S(i));      %S中每个输出数据流的完全传输时间
end
    
t=transt;
end

%*********************************************************
% II - check
function who=check(aa,bb,cc)
% 这个函数是为了计算Credit值，并且唯一目的是判断下一步该传输ABC中的谁。

%*************************
% 全局变量
global CreditA;
global CreditB;
global A;
global B;
global BE;

%*************************
% 功能部分

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

who=W;              %输出变量who将得到序列数据流的角标序号 
end

%*********************************************************
% III - trans
function tt=trans(who)      % 参数who被作为序列数据流的角标而调用

% trans只作为传输函数工作                            "trans" works only as a transmission function. 
% 它不仅仅选择传输时间，还可以计算credit的变化。      It not only selects the transmission time, but also calculates the change of credit. 
% 然后它就可以根据credit区分出被传输的对象。          Then it will distinguish the transmitted object according to Credit. 

%当有不止一个C级帧时：
    %首先传输C级帧。 
    %每次在C级帧传输完毕后，计算CreditA。如果CreditA>=0，则下一次传输A级帧；
    %A级帧每传输完一次，计算CA，如果CA>=0，则下一次传输A级帧；
                             %如果CA<0，计算CreditB：如果CB>=0，则下一次传输B级帧；
                                                   %如果CB<0，则传输C级帧。
    %B级帧每传输完一次，计算CA，如果CA>=0，则下一次传输A级帧；
                             %如果CA<0，计算CreditB：如果CB>=0，则下一次传输B级帧；
                                                   %如果CB<0，则传输C级帧。
                                                   
%*************************
% 全局变量
global CreditA;
global CreditB;
global aA_pos;
global aA_neg;
global aB_pos;
global aB_neg;
global C;
global pprr;

%*************************
% 功能部分
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