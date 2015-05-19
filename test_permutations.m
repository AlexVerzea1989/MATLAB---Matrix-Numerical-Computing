function r=test_permutations
n=20;
dim=[n,n];
reductions = {@QR,@LLLP,@SQRD,@VBLAST};
%reductions = {@LLLA,@LLLV,@VLLL};
nt=numel(reductions);
b=[0:0.05:0.7];
rept=1000;

rho=1;
r=zeros(nt,length(b));
increases=zeros(nt,length(b));

for k=1:length(b)
disp(b(k));
for i=1:rept
    RF=b(k);
    A=randILS2(dim,0.1,b(k),0,0);
    before_incr=0;
    for j=1:length(reductions)
    
            F=reductions{j};
            R=F(A,n);
            
            s=scmplx_beta2(R,rho);
            r(j,k)=r(j,k)+s;
            
            if before_incr==0
                before_incr=s;
            end
            if s>before_incr
                increases(j,k)=increases(j,k)+1;
            end
    end
end
end
r=r/rept;
names=cell(1,nt);
for i=1:nt
    names{i}=func2str(reductions{i});
end
plotS(r',names,b);
set(gca, 'yscale', 'log');

increases 

end

function R=QR(A,n)
[~,R]=qr(A,0);
end

function R=LLLP(A,n)
R=LLL_permute(A,n);
end

function R=SQRD(A,n)
R=QRP_PHH(A,n);
end

function R=VBLAST(A,n)
R=QRP_VBLAST(A,n);
end

function R=LLLA(A,n)
R=LLL(A,n);
end

function R=LLLV(A,n)
R=LLL(A,n);
R=QRP_VBLAST(R,n);
end

function R=VLLL(A,n)
R=QRP_VBLAST(A,n);
R=LLL(R,n);
end


% h=xlabel('$\delta$');
% set(h, 'Interpreter', 'latex');
% set(h,'FontSize',12);
% h=ylabel('Average $\bar{\eta}(\boldsymbol{R})$ over $100$ runs');
% set(h, 'Interpreter', 'latex');
% set(h,'FontSize',12);
% %ti = get(gca,'Position');
% %set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
% set(gca,'Position',[0.1049    0.1049    0.8623    0.8623]);
% p=get(gcf,'Position');
% p(3)=469;p(4)=290;
% set(gcf,'Position',p);
