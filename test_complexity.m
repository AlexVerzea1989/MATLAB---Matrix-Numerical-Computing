RF=0;
RH=0;
n=20;
dim=[n,n];
delta = [0.3:0.1:1];
b=[3.5,3,2.5,2]*10;
increase=zeros(length(b),length(delta));
before_incr=zeros(length(b),1);
r=zeros(length(delta),length(b));

for i=1:100
    fprintf('%d,',i);
    j=0;
    A=randILS2(dim,0.1,RH,RF,2);
    before_incr=zeros(length(b),1);
    for j=1:length(delta)
    
        if delta(j)<0.25
            [~,R]=qr(A);
        else
            R=LLL(A,n,delta(j));
        end
        
        for k=1:length(b)
            s=scmplx_beta2(R,b(k));
            r(j,k)=r(j,k)+s;
            
            if(s>before_incr(k))
                increase(k,j)=increase(k,j)+1;
            end
            before_incr(k)=s;
        end
        
        
    end
end
r=r/100;

plotS(r,{'$\rho=35$','$\rho=30$','$\rho=25$','$\rho=20$'},delta);
%plotS(r,{'$\rho=3.5$','$\rho=3.0$','$\rho=2.5$','$\rho=2.0$'},delta);
set(gca, 'yscale', 'log');
h=xlabel('$\delta$');
set(h, 'Interpreter', 'latex');
set(h,'FontSize',12);
h=ylabel('Average $\bar{\eta}(\boldsymbol{R})$ over $100$ runs');
set(h, 'Interpreter', 'latex');
set(h,'FontSize',12);
%ti = get(gca,'Position');
%set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
set(gca,'Position',[0.1049    0.1049    0.8623    0.8623]);
p=get(gcf,'Position');
p(3)=469;p(4)=290;
set(gcf,'Position',p);

increase'
