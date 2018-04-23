%% ��������Sequence
function y=Sequence(S)    % �����������  Ordered Output Function. 
% ����Ӧ����S1/S2/S3
% Seq����һ���������������������������ĽǱꡣ
% y����������е�˳�� 

%********************************************************
% ȫ�ֱ���
%global S1;          %S1=[1 0 0; 2 3 6; 0 0 0];
%global S2;          %S2=[1 2 3; 4 0 0; 5 0 0];
%global S3;          %S3=[2 0 0; 3 0 0; 6 0 0];%�����S��ʽ��Ҫ���ģ�������ܻ����ά�Ȳ���ȷ
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

%% �Ӻ���
%*********************************************************
% I - OUT
function y=OUT(S)    
% S��һ������˿ڣ�SΪS1/S2/S3... ȫ�ֱ���I�ǽǱ꣬�������ȼ���S is an input port, S=S1/S2/S3... I is subscript, not priority!
% Try to unify the function of Sequence into one module, and run it three times with different parameters.

%*************************
% ȫ�ֱ���
global I;
global Pr;

%*************************
%�ڽ���������˿ڴ��������������Ϊ�������롣 
pr=0;
A=size(S);                  %�õ�S�����룩���С�����Ŀ�����д���������˿ڣ��д���ĳһ���˿�����������������������
k=1;
for i=1:A(1)                %A(1)��������˿�����
    for j=1:A(2)            %A(2)������༸������������
        Seqm_n(k)=S(i,j);   %������ԭʼ��������Sת��Ϊ����������Seqm_n����Ϊһ�����������ۼ�����ڼ������ڣ�һ�ζ�ֻ�ܴ���һ�����ݡ���
        k=k+1;
    end
end

%��Seqm_n�еķ�0Ԫ����ȡ�����õ�������������������������Seq��Seq�е�����Ϊ�Ǳ꣩
k=1;
for i=1:length(Seqm_n)
    if (Seqm_n(i)~=0)
        temp=Seqm_n(i);     
        Seq(k)=temp;
        k=k+1;
    end
end

%������������������������Seq��Ԫ�ص����ȼ���ȡ�����õ�����pr
for i=1:length(Seq)
    pr(i)=Pr(Seq(i));     
end

%��������������ͨ��liner_i���õ�������������У��ѿ�������ӳٶ��󣩵ĽǱ�����ȼ����������Ϊs������
s=liner_i(pr,Seq,I);        %ͨ������������������о�����         % s��һ��������s=[pr;Seq]
SIZES_m=size(s);            % SIZES_m����s���С�����Ŀ����s�����̶�Ϊ2������Ϊ�����������������

for i=1:SIZES_m(2)          %SIZES_m(2)����Seq����������s��������Ҳ�����������  
    S_mn(i)=s(2,i);         %S_mn�ǿ�������ӳٶ����������У�������Ϊ������еĽǱ꣨���ǰ�s�ĵ�2����ȡ�����ˣ�
end

%for i=1:SIZESm(2)          % ��s����ȡ���ȼ����У�����1��
%    Pr_mn(i)=s(1,i);
%end

%*************************
% ���
y=S_mn;                     %yΪ������нǱ�
end

%*********************************************************
% II - liner_i
function y=liner_i(pr,Seq,I)
% ����������ӳٶ�������ȼ���������ͬʱ���������ĽǱ�����
% I�ǽǱ���������ȼ�

%*************************
% ȫ�ֱ���
global Pr;

%*************************
%������BE�������ŵ��������е���ǰ��
for i=1:length(pr)          %�ӵ�һ֡�����һ֡ 
    m=1;
    if pr(i)==1             %��������κ�BE�� If there exists any BE flow, 
        temp=pr(i);         %�ͽ���һ֡�ʹ�ͷ���ο�ʼ����һ֡���н������ȸ���1֡�����ٸ���2֡�����ٸ���3֡������
        pr(i)=pr(m);
        pr(m)=temp;
        temp1=Seq(i);
        Seq(i)=Seq(m);
        Seq(m)=temp1;
        m=m+1;
    end
end                         %����BE��������������������ǰ��

%*************************
%��Class B���ŵ������������
for i=1:length(pr)          %�ӵ�һ֡�����һ֡
    k=length(pr);
    if pr(i)==2             %�����B�� 
        temp=pr(i);         %�ͽ�����һ֡�ʹ�������ε�����ǰ�ߵ���һ֡
        pr(i)=pr(k);
        pr(k)=temp;
        temp1=Seq(i);
        Seq(i)=Seq(k);
        Seq(k)=temp1;
        k=k-1;
    end
end                         %����B���������������������

%*************************
% Put the Researched Data Stream at the Last Part of the Class Which Is Being Studied. ��ĳһ�������ȼ��ڣ������о��������ŵ����
a=0;
b=0;
c=0;
for i=1:length(pr)          % Collect the number of each kind of priorities of output data stream. �ռ�����������и����ȼ��ĸ���
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