function fs=testReduct_box(rk,rpt,interval,noise,type)
%global sigma;
global cpu_time_0;
global Rho_2;
global delta2
Rho_2=inf;
global ILS_RSD;
ILS_RSD=nan;
%tests={@testPLLL_BC_IGT,@testBC,@testCCPX};
tests={@testISET_IGT2,@testBC};
%tests={@testBC,@testNBC};
%tests={@testBC,@testCCP};
%tests={@testSQRD,@testBC};
%tests={@testNew3,@testCCP};
%objs={@objFastRsd,@objBabaiRsd,@objRoundRsd,@objILSRsd,@objTotalTime;@objSearchTime};objReductionTime
%objs={@objISF2,@objISF2PLLL};
objs={@objSearchTime};

nt=numel(tests);
no=numel(objs);
nnnsss=max([size(rk,2),length(noise),length(interval)]);
redtime=zeros(nt,1);
xnoise=false;
xinterval=false;
if(length(noise)==1)
    noise=noise*ones(1,nnnsss);
end
if(size(rk,2)==1)
    rk=repmat(rk,1,nnnsss);
    xnoise=true;
end
if(size(interval,2)==1)
    interval=interval*ones(1,nnnsss);
else
    xinterval=true;
end
fs=zeros(nnnsss,nt*no);
%snr=zeros(size(rk,2),1);
%M=2^(2*interval);
%disp([num2str(M),'-QAM']);
%disp(['SNR: ',num2str(10*log10(((M-1)/6)/noise^2))]);
%sigma = sqrt(noise);


for n=1:nnnsss
    rn=rk(:,n)';
    %disp(rn)
    for j=1:rpt
        
        [A,y,l,u,~,arn]=randILS_box(type,rn,noise(n),interval(n));
        fprintf(['\n',num2str(arn),'/',num2str(noise(n)),'/',num2str(interval(n)),':',num2str(j),'-']);
        %snr(n)=10*log10(((M-1)/6)/noise^2);
        %disp(cond(A));
        delta2=noise(n)^2;
        for i=1:nt
            testi=tests{i};
            fprintf(func2str(testi));
            fprintf('(');
            Rho_2=inf;
            tic;
            [Ay2,Z2,l2,u2]=testi(A,y,l,u,arn(2));
            cpu_time_0=toc;
            fprintf([num2str(cpu_time_0),'|']);
            redtime(i)=redtime(i)+cpu_time_0;
            %disp(testi)
            %disp(diag(Ay2)')
            for k=1:no
                obji=objs{k};
%                 if isequal(tests{i},@testVBLAST) && noise>=0.5 && arn(2)>=9
%                     f=nan;
%                 else
                    f=obji(Ay2,Z2,l2,u2,arn);
%                 end
                %                if(f>2e-3 && i==2)
                %                    keyboard;
                %                end
                fs(n,(k-1)*nt+i)=fs(n,(k-1)*nt+i)+f;
                fprintf([num2str(f),',']);
            end
            fprintf(')-');
        end
    end
end
fs=fs/rpt;
%[xaxis,ix]=sort(rk);
names=cell(1,nt*no);
for i=1:nt
    name=func2str(tests{i});
    name=name(5:length(name));
    for j=1:no
        objName=func2str(objs{j});
        names{nt*(j-1)+i}=[name,':',objName(4:length(objName))];
    end
end

if no>1
    if(xinterval)
        xlabel('interval:2^n');
        plotS(fs,names,interval,nt);
    elseif(xnoise)
        xlabel('noise');
        plotS(fs,names,noise,nt);
    else
        xlabel('dimension');
        plotS(fs,names,2*rk(1,:),nt);
    end
end

for obji=1:no
    fsi=fs(:,((obji-1)*nt+1):(obji*nt));
    namei={names{((obji-1)*nt+1):(obji*nt)}};
figure;hold on;ylabel(['Average over ',num2str(rpt),' runs']);
grid on;
if(xinterval)
    xlabel('interval:2^n');
    plotS(fsi,namei,interval,nt);
elseif(xnoise)
    xlabel('noise');
    plotS(fsi,namei,noise,nt);
else
    xlabel('dimension');
    plotS(fsi,namei,2*rk(1,:),nt);
end
%set(gca,'YScale','log');
disp(redtime)
end
end

function [Ay,Z,l,u]=testSQRD(A,y,l,u,n)
[Ay,Z,l,u] = QRZ_sqrd([A,y],l,u,n);
end

function [Ay,Z,l,u]=testVBLAST(A,y,l,u,n)
[A,y,l,u,~,Z]=QRZ_vblast(A,y,l,u);
Ay=[A,y];
end

function [Ay,Z,l,u]=testNew4(A,y,l,u,n)
[Ay,Z,l,u]=QRZ_CCP4([A,y],l,u,n);
end

function [Ay,Z,l,u]=testNew3(A,y,l,u,n)
[Ay,Z,l,u]=QRZ_CCP4([A,y],l,u,n);
end

function [Ay,Z,l,u]=testCCP(A,y,l,u,n)
[Ay,Z,l,u]=QRZ_CCP([A,y],l,u,n);
end

function [Ay,Z,l,u]=testNone(A,y,l,u,n)
Z=eye(n);
Ay=QRZ_HH([A,y],n);
end

function [Ay,Z,l,u]=testXPLLL2(A,y,l,u,n)
[Ay,Z,l,u]=XPLLL2([A,y],l,u,n);
end

function [Ay,Z,l,u]=testXPLLL(A,y,l,u,n)
[Ay,Z,l,u]=XPLLL([A,y],l,u,n);
end


function [Ay,Z,l,u]=testNoIGT(A,y,l,u,n)
[Ay,Z,l,u]=PLLL3([A,y],l,u,n,zeros(n,1));
end
%
%
% function [Ay,Z,l,u]=testNew3(A,y,l,u,n)
% [Ay,Z,l,u]=QRZ_CX_box([A,y],l,u,n,n+1,3);
% end
%
% function [Ay,Z,l,u]=testNew4(A,y,l,u,n)
% [Ay,Z,l,u]=QRZ_CX_box([A,y],l,u,n,n+1,4);
% end

function [Ay,Z,l,u]=testRepeatNew(A,y,l,u,n)
Ay=[A,y];
Z=eye(n);
i=10;
while i>0
    [Ay,Z2,l,u]=QRZ_HH_box(Ay,l,u,n);
    if sum(diag(Z2))==n
        break
    end
    Z=Z2*Z;
    i=i-1;
end
end

function [Ay,Z,l,u]=testCH(A,y,l,u,n)
[R,Z,y,l,u] = QRZ_CH(A,y,l,u);
Ay=[R,y];
end

%CH
function [Ay,Z,l,u]=testBC(A,y,l,u,n)
%[A,Z,y,l,u]=QRZ_BC(A,y,l,u);
[A,y,Z,l,u]=QRP_BC(A,y,l,u);
Ay=[A,y];
end

function [Ay,Z,l,u]=testBC_ISET(A,y,l,u,n)
[A,y,p,l,u]=QRZ_BC_ISET(A,y,l,u);
Ay=[A,y];
Z=eye(n);
Z=Z(p,:);
end


%CH
function [Ay,Z,l,u]=testNBC(A,y,l,u,n)
[A,y,Z,l,u]=QRP_NBC(A,y,l,u);
Ay=[A,y];
end

function [Ay,Z,l,u]=testCCPX(A,y,l,u,n)
[A,y,Z,l,u] = QRZ_CCPX(A,y,l,u);
Ay=[A,y];
end

%LLL+(CH+IGT)
function [Ay,Z,l,u]=testPLLL_BC_IGT(A,y,l,u,n)
[A,y,Z,l,u] = QRZ_XIE(A,y,l,u,0);
Ay=[A,y];
end

%LLL+(CH+IGT)
function [Ay,Z,l,u]=testPLLL_BC_IGT2(A,y,l,u,n)
%[A,y,Z,l,u] = QRZ_XIE(A,y,l,u,-2);
[A,y,~,l,u] = QRZ_XIE(A,y,l,u,0);
[A,y,Z,l,u] = QRZ_XIE(A,y,l,u,0);
Ay=[A,y];
end

%New method
function [Ay,Z,l,u]=testNew(A,y,l,u,n)
[A,y,Z,l,u] = QRZ_XIE(A,y,l,u,4);
Ay=[A,y];
end

%CH+IGT
function [Ay,Z,l,u]=testBC_IGT(A,y,l,u,n)
[A,y,Z,l,u] = QRZ_XIE(A,y,l,u,1);
Ay=[A,y];
end

%LLL+New method
function [Ay,Z,l,u]=testPLLL_New(A,y,l,u,n)
[A,y,Z,l,u] = QRZ_XIE(A,y,l,u,3);
Ay=[A,y];
end

%LLL+CH method
function [Ay,Z,l,u]=testPLLL_BC(A,y,l,u,n)
[A,y,Z,l,u] = QRZ_XIE(A,y,l,u,2);
Ay=[A,y];
end

%LLL method
function [Ay,Z,l,u]=testPLLL(A,y,l,u,n)
[Ay,Z] = PLLL([A,y],n);
end

%UPDATE method
function [Ay,Z,l,u]=testUPDATE(A,y,l,u,n)
%[A,y,l,u,p]=QRZ_vblast(A,y,l,u);
% [Q,A]=qr(A);
% y=Q'*y;
[A,y,l,u,p]=QRP_NBC(A,y,l,u);
inactive=zeros(size(l));
[R,y,l,u,Z]=QRP_UPDATE(A,y,l,u,inactive);
Z=Z(p,:);
Ay=[R,y];
end

%LLL method
function [Ay,Z,l,u]=testISET_IGT(A,y,l,u,n)
[A,y,Z,l,u]=QRZ_ISET_IGT(A,y,l,u,false,true,2);
Ay=[A,y];
end

%LLL method
function [Ay,Z,l,u]=testISET_IGT2(A,y,l,u,n)
[A,y,Z,l,u]=QRZ_ISET_IGT(A,y,l,u,false,true,13);
Ay=[A,y];
end

%LLL method
function [Ay,Z,l,u]=testISET_IGT3(A,y,l,u,n)
[A,y,Z,l,u]=QRZ_ISET_IGT(A,y,l,u,false,true,1);
Ay=[A,y];
end

%LLL method
function [Ay,Z,l,u]=testPLLLBC(A,y,l,u,n)
[Ay,Z] = PLLL([A,y],n);
[A,y,l,u,p]=QRP_NBC(Ay(1:n,1:n),Ay(1:n,n+1),l,u);
Z=Z(p,:);
Ay=[A,y];
end

function [Ay,Z,l,u]=testTSD(A,y,l,u,n)

Proj_y=y+A*ones(n,1);
Proj_A=2*A;

Proj_d=-(Proj_y'*Proj_A)';
Proj_x=gradproj(0.5*ones(n,1),Proj_A,Proj_d,1,0);
Proj_x_2=2*round(Proj_x)-1;

Proj_gamma=norm(y-A*Proj_x_2);

[Q,R,Z]=Reduction(A,y,Proj_gamma);
%A*P=Q*R
y=Q'*y;
l=round(Z\l);
u=round(Z\u);
Ay=[R,y];
end

function f=objSearchComplexity(Ay,Z,l,u)
search_CSD_box(Ay,l,u);
f=set_count(0,2);
end

function f=objScmplx(Ay,Z,l,u)
f=scmplx_beta(Ay,10);
end

function f=objTotalTime(Ay,Z,l,u,dim)
global cpu_time_0;
f=objSearchTime(Ay,Z,l,u,dim);
f=f+cpu_time_0;
end

function f=objReductionTime(Ay,Z,l,u,dim)
global cpu_time_0;
f=cpu_time_0;
end

function f=objSearchTime(Ay,Z,l,u,dim)
global Rho_2;
global ILS_RSD;
m=dim(1);
n=dim(2);
if m<n
    R=Ay(1:m,1:n);
    y=Ay(1:m,n+1);
    tic
    [~,ILS_RSD]=search_UBILS(R,y,l,u,1);
    f=toc;
else
    tic
    [~,ILS_RSD]=search_CSD_box(Ay,l,u,Rho_2);
    f=toc;
end
%Rho_2=inf;
end

function f=objRho(Ay,Z,l,u,dim)
global Rho_2;
f=Rho_2;
%Rho_2=inf;
end

function f=objNumISet(Ay,Z,l,u,dim)
f=sum(isinf(l)&isinf(u));
end

function f=objISF(Ay,Z,l,u,dim)
n=dim(1);
R=Ay(1:n,1:n);
y=Ay(:,n+1);
f = sum(ISF(R,y,l,u));
end

function f=objISFT(Ay,Z,l,u,dim)
n=dim(1);
R=Ay(1:n,1:n);
y=Ay(:,n+1);
tic
ISF(R,y,l,u);
f=toc;
end

function f=objISF1(Ay,Z,l,u,dim)
n=dim(1);
R=Ay(1:n,1:n);
y=Ay(:,n+1);
f = sum(ISF1(R,y,l,u));
end

function f=objISF1T(Ay,Z,l,u,dim)
n=dim(1);
R=Ay(1:n,1:n);
y=Ay(:,n+1);
tic
ISF1(R,y,l,u);
f=toc;
end

function f=objISF2(Ay,Z,l,u,dim)
n=dim(1);
R=Ay(1:n,1:n);
y=Ay(:,n+1);
f = sum(ISF2(R,y,l,u,1));
end

function f=objISF2w(Ay,Z,l,u,dim)
n=dim(1);
R=Ay(1:n,1:n);
y=Ay(:,n+1);
f = sum(ISF2(R,y,l,u,0));
end

function f=objISF2PLLL(Ay,Z,l,u,dim)
n=dim(1);
R=Ay(1:n,1:n);
y=Ay(:,n+1);
r2=PLLL_residual(Ay,l,u);
f = sum(ISF2(R,y,l,u,r2));
end

function f=objISF2T(Ay,Z,l,u,dim)
n=dim(1);
R=Ay(1:n,1:n);
y=Ay(:,n+1);
tic
ISF2(R,y,l,u);
f=toc;
end

function f=objBabaiRsd(Ay,Z,l,u,dim)
n=dim(2);
%[~,f] = search_babai_box(Ay(n-1:end,n-1:end),l(n-1:end),u(n-1:end));
[~,f] = search_babai_box(Ay,l,u);
end

function f=objINUM(Ay,Z,l,u,dim)
f=sum((u-l)==0);
end

function f=objRoundRsd(Ay,Z,l,u,dim)
x=Ay(:,1:dim(2))\Ay(:,dim(2)+1);
x=max(l,min(u,round(x)));
r=Ay*[x;-1];
f=r'*r;
end

function f=objILSRsd(Ay,Z,l,u,dim)
global ILS_RSD;
if isnan(ILS_RSD);
[~,f] = search_CSD_box(Ay,l,u);
else
    f=ILS_RSD;
end
ILS_RSD=nan;
%f2=objBabaiRsd(Ay,Z,l,u,dim);
% if f2>2*f
%     keyboard;
% end
end

function f=objSearchTimeUnbounded(Ay,Z,l,u,dim)
tic
search_CSD(Ay);
f=toc;
end
