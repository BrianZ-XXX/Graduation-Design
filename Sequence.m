%% 主函数：Sequence
function y=Sequence(S)    % 有序输出函数  Ordered Output Function. 
% 输入应该是S1/S2/S3
% Seq代表一个向量，并且其内容是数据流的角标。
% y代表输出序列的顺序 

%********************************************************
% 全局变量
%global S1;          %S1=[1 0 0; 2 3 6; 0 0 0];
%global S2;          %S2=[1 2 3; 4 0 0; 5 0 0];
%global S3;          %S3=[2 0 0; 3 0 0; 6 0 0];%下面的S格式需要更改，否则可能会出现维度不正确
global Pr;
global SS
%global S11;
%global S21;
%global Pr31;

SS=OUT(S);
%S21=OUT(S2);
%S31=OUT(S3);
for i=1:length(SS)
    pr(i)=Pr(SS(i));
end
y=[SS;pr];
end

%% 子函数
%*********************************************************
% I - OUT
function y=OUT(S)    
% S是一个输入端口，S为S1/S2/S3... 全局变量I是角标，不是优先级！S is an input port, S=S1/S2/S3... I is subscript, not priority!
% Try to unify the function of Sequence into one module, and run it three times with different parameters.

%*************************
% 全局变量
global I;
global Pr;

%*************************
%在交换机输入端口处将并行输入调整为串行输入。 
pr=0;
A=size(S);                  %得到S（输入）的行、列数目。（行代表几个输入端口，列代表某一个端口中输入最多的数据流数量）
k=1;
for i=1:A(1)                %A(1)代表输入端口数量
    for j=1:A(2)            %A(2)代表最多几个数据流数量
        Seqm_n(k)=S(i,j);   %将无序原始键入数据S转变为无序串行输入Seqm_n（因为一个交换机无论几个入口几个出口，一次都只能传递一个数据。）
        k=k+1;
    end
end

%将Seqm_n中的非0元素提取出，得到真正的无序串行输入数据向量Seq（Seq中的内容为角标）
k=1;
for i=1:length(Seqm_n)
    if (Seqm_n(i)~=0)
        temp=Seqm_n(i);     
        Seq(k)=temp;
        k=k+1;
    end
end

%将真正的无序串行输入数据流Seq内元素的优先级提取出来得到向量pr
for i=1:length(Seq)
    pr(i)=Pr(Seq(i));     
end

%将串行输入数据通过liner_i并得到有序串行输出序列（已考虑了最坏延迟对象）的角标和优先级，输出数量为s的列数
s=liner_i(pr,Seq,I);        %通过第三个参数决定最坏研究对象         % s是一个向量，s=[pr;Seq]
SIZES_m=size(s);            % SIZES_m代表s的行、列数目。（s行数固定为2，列数为输出数据流的数量）

for i=1:SIZES_m(2)          %SIZES_m(2)代表Seq的数量，即s的列数，也即输出数量。  
    S_mn(i)=s(2,i);         %S_mn是考虑了最坏延迟对象的输出序列，其内容为输出序列的角标（就是把s的第2行提取出来了）
end

%for i=1:SIZESm(2)          % 从s中提取优先级序列，即第1行
%    Pr_mn(i)=s(1,i);
%end

%*************************
% 输出
y=S_mn;                     %y为输出序列角标
end

%*********************************************************
% II - liner_i
function y=liner_i(pr,Seq,I)
% 将考虑了最坏延迟对象的优先级序列排序，同时将数据流的角标排序。
% I是角标而不是优先级

%*************************
% 全局变量
global Pr;

%*************************
%将所有BE数据流放到整个序列的最前面
for i=1:length(pr)          %从第一帧到最后一帧 
    m=1;
    if pr(i)==1             %如果存在任何BE流 If there exists any BE flow, 
        temp=pr(i);         %就将这一帧和从头依次开始的那一帧进行交换（先跟第1帧换，再跟第2帧换，再跟第3帧换…）
        pr(i)=pr(m);
        pr(m)=temp;
        temp1=Seq(i);
        Seq(i)=Seq(m);
        Seq(m)=temp1;
        m=m+1;
    end
end                         %所有BE流被换到了整个序列最前方

%*************************
%将Class B流放到整个序列最后
for i=1:length(pr)          %从第一帧到最后一帧
    k=length(pr);
    if pr(i)==2             %如果有B流 
        temp=pr(i);         %就交换这一帧和从最后依次倒着往前走的那一帧
        pr(i)=pr(k);
        pr(k)=temp;
        temp1=Seq(i);
        Seq(i)=Seq(k);
        Seq(k)=temp1;
        k=k-1;
    end
end                         %所有B流被换到了整个序列最后方

%*************************
% Put the Researched Data Stream at the Last Part of the Class Which Is Being Studied. 在某一大类优先级内，将被研究数据流放到最后
a=0;
b=0;
c=0;
for i=1:length(pr)          % Collect the number of each kind of priorities of output data stream. 收集输出数据流中各优先级的个数
    if pr(i)==3
        a=a+1;
    elseif pr(i)==2
        b=b+1;
    else
        c=c+1;
    end
end

k=c+a;
if Pr(I)==3                 % If the studied object is Class A.
    for i=(c+1):(c+a)
        if Seq(i)==I        % Meeting the worst delay object which is studied. 
            temp=pr(i);     % Exchange this frame and the last frame in this class. 
            pr(i)=pr(k);
            pr(k)=temp;
            temp1=Seq(i);
            Seq(i)=Seq(k);
            Seq(k)=temp1;
        end
    end
end

k=length(pr);
if Pr(I)==2                 % If the studied object is Class B.
    for i=(a+c+1):length(pr)
        if Seq(i)==I        % Meeting the worst delay object which is studied. 
            temp=pr(i);     % Exchange this frame and the last frame in this class.
            pr(i)=pr(k);
            pr(k)=temp;
            temp1=Seq(i);
            Seq(i)=Seq(k);
            Seq(k)=temp1;
        end
    end
end

%*************************
% the Non-researched Data of SR Stream Will Be at the Middle Section (Autocomplete)

%*************************
% Output
y=[pr;Seq];
end