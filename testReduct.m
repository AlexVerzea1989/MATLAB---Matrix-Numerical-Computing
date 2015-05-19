function [fs,xaxis]=testReduct(rk,rpt,econd,noise)
global sigma;
global cpu_time;
tests={@testPLLL2,@testPLLL};
%objs={@objFastRsd,@objBabaiRsd,@objILSRsd};
%objs={@objSearchTime,@objRsd,@objSearchTimeFast3,@objRsd,@objSearchTimeFast2,@objRsd};
%objs={@objSearchTime,@objSearchTimeFast,@objSearchTimeFast2};
objs={@objToc};
nt=numel(tests);
no=numel(objs);
fs=zeros(length(rk),nt*no);
sigma = sqrt(noise);
for n=1:length(rk)
    rn=rk(n);
    disp(rn)
    for j=1:rpt
        [A,y]=randILS(rn,rn,econd,noise);
        rn=size(A,2);
        for i=1:nt
            testi=tests{i};
            tic;
            [Ry,Z]=testi(A,y,rn);
            %disp(Z)
            cpu_time=toc;
            for k=1:no
                obji=objs{k};
                f=obji(Ry,rn);
                fs(n,(k-1)*nt+i)=fs(n,(k-1)*nt+i)+f;
            end
        end
    end
end
fs=fs/rpt;

[xaxis,ix]=sort(rk);
names=cell(1,nt*no);
for i=1:nt
    name=func2str(tests{i});
    name=name(5:length(name));
    for j=1:no
        objName=func2str(objs{j});
        names{nt*(j-1)+i}=[name,':',objName(4:length(objName))];
    end
end
fs=fs(ix,:);
figure;hold on;ylabel(['Average over ',num2str(rpt),'runs']);xlabel('Dimension');
grid on;plotS(fs,names,xaxis,nt);

end

function [Ry,Z]=testNew(R,y,n)
[Ry,Z] = XHReduce([R,y],n);
end

function [Ry,Z]=testNewSW(R,y,n)
[Ry,Z]=testSW(R,y,n);
[Ry,Z2] = XHReduce(Ry,n);
Z=Z*Z2;
end

function [Ry,Z]=testNewVB(R,y,n)
[Ry,Z]=testVB(R,y,n);
[Ry,Z2] = XHReduce(Ry,n);
Z=Z*Z2;
end

function [Ry,Z]=testVB(R,y,n)
[Ry,Z]=QRZ_PHH([R,y],n);
end

function [Ry,Z]=testSW(R,y,n)
[R,y,Z] = otherConstrainedReduction(R,y);
Ry=[R,y];
end

function [Ry,Z]=testNone(R,y,n)
[Ry,Z] = QRZ_HH([R,y],n);
end

function [Ry,Z]=testPLLL(R,y,n)
[Ry,Z] = PLLL([R,y],n);
end

function [Ry,Z]=testPLLL2(R,y,n)
[Ry,Z] = PLLL2([R,y],n);
end

function f=objFastRsd(Ry,n)
[~,f] = search_fast(Ry,n);
end

function f=objBabaiRsd(Ry,n)
[~,f] = search_babai(Ry,n);
end

function f=objILSRsd(Ry,n)
[~,f] = search_CSD(Ry,n);
end

function f=objToc(Ry,n)
global cpu_time;
f=cpu_time;
end

function f=objComplexity1(Ry,n)
f=set_count(0);
end

function f=objComplexity2(Ry,n)
f=set_count(0,2);
end

function f=objComplexityS1(Ry,n)
set_count(0);
set_count(0,2);
[Z,rsd]=search_CSD(Ry,inf);
f=set_count(0,2);
end

function f=objPLLLReduced(Ry,n)
f=isPLLLReduced(Ry,1,n);
end

function f=objComplexityS2(Ry,n)
set_count(0);
set_count(0,2);
[Z,rsd]=search_CSD(Ry,-inf);
f=set_count(0,2);
end

function f=objComplexityT(Ry,n)
[z,beta,f]=scmplx(Ry,n);
end

function f=objSRate(Ry,n)
global sigma;
f=srate_babai(Ry,sigma,n);
end

function f=objSearchTime(Ry,n)
global Z;
tic;
[Z,rsd]=search_CSD(Ry);
f=toc;
%disp('normal:');
%disp(Z');
end

function f=objSearchTimeFast(Ry,n)
global Z;
tic;
[Z,rsd]=search_fast(Ry,n,n+1,n);%n/3
[Z,rsd]=search_CSD(Ry,rsd*1.0001);
f=toc;
end

function f=objRsd(Ry,n)
global Z;
f= sum((Ry*[Z;-1]).^2);
end

function f=objSearchTimeFast2(Ry,n)
global Z;
tic;
[Z,rsd]=search_fast_n(Ry,n,n+1,n);%n/3
[Z,rsd]=search_CSD(Ry,rsd*1.0001);
f=toc;
end
% 
% function f=objError(Ry,n,A,y,Z)
%             [~,R]=qr(A*Z);
%             D=R-Ry(1:rn,1:rn);
%             norm(D/Z)
% end

function f=objSearchTimeFast3(Ry,n)
global Z;
tic;
[Z,rsd]=search_fast2(Ry,n,n+1,n);%n/3
f=toc;
end
