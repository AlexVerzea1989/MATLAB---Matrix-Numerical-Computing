function test(rk,rpt,econd,id)
%tests={@testLLL,@testELLL,@testPLLL};
tests={@testPLLL,@testPLLL2};
time=zeros(length(rk),length(tests));
%flops=zeros(length(rk),length(tests));
errs=zeros(length(rk),length(tests));
% stime=zeros(length(rk),length(tests));
% btime=zeros(length(rk),length(tests));
biterr=zeros(length(rk),length(tests));
%biterrsd=zeros(length(rk),length(tests));
nt=numel(tests);
delta=4/3;
for n=1:length(rk)
    rn=rk(n);
    disp(rn)
    for j=1:rpt
%         disp(j)
        [A,y,x]=randILS(rn,rn,econd);
        Ay=[A,y];
        AI=[A,eye(rn)];
        for i=1:nt
            testi=tests{i};
            %-------time----
            tic;
            [Ry,Z]=testi(Ay,rn,delta);
            R=Ry(:,1:rn);
            y=Ry(:,rn+1);
            t=toc;
            time(n,i)=time(n,i)+t;
%            flops(n,i)=flops(n,i)+get_count();
            %-------err-------
             [RQ]=testi(AI,rn,delta);
             Q=RQ(:,rn+1:2*rn);
             err=norm(A-Q'*R/Z)/norm(A);
             errs(n,i)=errs(n,i)+err;
            %-------SD-----
%             tic;
%             [zs,rs]=search_NSD(R,y);
% %             t=toc;
%             sx=Z*zs;
%             biterrsd(n,i)=biterrsd(n,i)+sum((sx-x)~=0);
%             stime(n,i)=stime(n,i)+t;
            %-------Babai---
%             tic;
            [zb,rb]=search_babai(Ry,rn);
%             t=toc;
            bx=Z*zb;
%             if norm(zs-zb,inf)>0.000001
%                 disp([zs,zb]);
%                 disp([rs,rb]);
%             end
%             btime(n,i)=btime(n,i)+t;
            biterr(n,i)=biterr(n,i)+sum((bx-x)~=0);
        end
    end
    time(n,:)=time(n,:)/rpt;
    errs(n,:)=errs(n,:)/rpt;
%     stime(n,:)=stime(n,:)/rpt;
%     btime(n,:)=btime(n,:)/rpt;
    %flops(n,:)=flops(n,:)/rpt;
    biterr(n,:)=biterr(n,:)/rpt/rn;
    %biterrsd(n,:)=biterrsd(n,:)/rpt/rn;
end
[xaxis,ix]=sort(rk);
names=cell(size(tests));
for i=1:numel(tests)
    names{i}=func2str(tests{i});
    names{i}=names{i}(5:length(names{i}));
end
save(['result/test',id,'.mat']);
% plotS(time(ix,:),names,xaxis);ylabel('reduction
% time');xlabel('Dimension');
% plotS(biterrsd(ix,:),names,xaxis);ylabel('Average bit error rate(SD)');xlabel('Dimension');grid on;saveas(gcf,['result/sberr',id,'.fig']);
% plotS(stime(ix,:),names,xaxis);ylabel('SD time');xlabel('Dimension');
% plotS(btime(ix,:),names,xaxis);ylabel('Babai time');xlabel('Dimension');
% 
gsize='-S';
%figure;hold on;ylabel('Average number of flops');xlabel('Dimension');grid on;plotS(flops(ix,:),names,xaxis,1,['result/flops',id],gsize);
figure;hold on;ylabel('Average running time');xlabel('Dimension');grid on;plotS(time(ix,:),names,xaxis,1,['result/time',id],gsize);
figure;hold on;ylabel('Average symbol error rate(Babai)');xlabel('Dimension');grid on;plotS(biterr(ix,:),names,xaxis,1,['result/bberr',id],gsize);
%figure;hold on;ylabel('Average symbol error rate(SD)');xlabel('Dimension');grid on;plotS(biterrsd(ix,:),names,xaxis,1,['result/bberrsd',id],gsize);
figure;hold on;ylabel('Backward errors');set(gca,'Yscale', 'log');grid on;xlabel('Dimension');plotS(errs(ix,:),names,xaxis,1,['result/err',id],gsize);
% gsize='-M';
% figure;hold on;ylabel('Average number of flops');xlabel('Dimension');grid on;plotS(flops(ix,:),names,xaxis,['result/flops',id],gsize);
% figure;hold on;ylabel('Average bit error rate(Babai)');xlabel('Dimension');grid on;plotS(biterr(ix,:),names,xaxis,['result/bberr',id],gsize);
% figure;hold on;ylabel('Backward errors');set(gca,'Yscale', 'log');grid on;xlabel('Dimension');plotS(errs(ix,:),names,xaxis,['result/err',id],gsize);
end

function [A,Z,invZ]=testLLL(A,n,delta)
    [A,Z,invZ]=LLL(A,n,delta);
end

function [A,Z,invZ]=testELLL(A,n,delta)
    [A,Z,invZ]=ELLL(A,n,delta);
end

function [A,Z,invZ]=testPLLL_1(A,n,delta)
    %stable,permute
    [A,Z,invZ]=PLLL(A,n,delta,@QRZ_PHH,@permu_Givens,1);
end

function [A,Z,invZ]=testPLLL_A(A,n,delta)
    %stable
    [A,Z,invZ]=PLLL(A,n,delta,@QRZ_GS,@permu_GS);
end

function [A,Z,invZ]=testPLLL_B(A,n,delta)
    %permute
    [A,Z,invZ]=PLLL(A,n,delta,@QRZ_HH,@permu_Givens);
end

function [A,Z]=testPLLL(A,n,delta)
    %permute
    [A,Z]=PLLL(A,n,delta,@QRZ_PHH);
end

function [A,Z]=testPLLL2(A,n,delta)
    %permute
    [A,Z]=PLLL2(A,n,delta,@QRZ_PHH);
end

function [A,Z,invZ]=testPLLL_N(A,n,delta)
    %non
    [A,Z,invZ]=ELLL(A,n,delta,@QRZ_HH,@permu_Givens);
end
