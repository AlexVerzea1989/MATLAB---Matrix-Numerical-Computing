function c=set_count(num,i)
global g_c;
if nargin<1
    num=0;
end
if nargin<2
    i=1;
end
if i<=numel(g_c)
    c=g_c(i);
else
    c=0;
end
g_c(i)=num;
