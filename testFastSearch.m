function testFastSearch(rk,rpt,type,noise)
t=zeros(3,length(rk));
rsd=zeros(3,length(rk));
errorbits=zeros(2,length(rk));
for n=1:length(rk)
    rn=rk(n);
    disp(rn)
    for j=1:rpt
%         disp(j)
        [A,y]=randAYX(type,rn,noise);
        Ay=[A,y];
        [Ry]=PLLL(Ay,rn);
        tic;
        [z1,rsd1]=search_fast_n(Ry);
        t1=toc;
        tic;
        [z2,rsd2]=search_CSD(Ry,1.0001*rsd1);
        t2=toc;
        tic;
        [z3,rsd3]=search_CSD(Ry);
        t3=toc;
        [z4,rsd4]=search_fast2(Ry);
        t(:,n)=t(:,n)+[t1;t2;t3];
        rsd(:,n)=rsd(:,n)+[rsd4;rsd1;rsd3];
        errorbits(:,n)=errorbits(:,n)+[sum(z1~=z3);sum(z2~=z3)];
    end
end
t=t/rpt;
t(2,:)=t(1,:)+t(2,:);
rsd=rsd/rpt;
errorbits(1,:)=errorbits(1,:)/rpt./rk;
errorbits(2,:)=errorbits(2,:)/rpt./rk;

figure;hold on;ylabel(['Average over ',num2str(rpt),' runs']);xlabel('Dimension');
grid on;title('search time');plotS(t',{'fast','fast+SD','SD'},rk);
figure;hold on;ylabel(['Average over ',num2str(rpt),' runs']);xlabel('Dimension');
grid on;title('residule');plotS(rsd',{'?','fast','SD'},rk);
% figure;hold on;ylabel(['Average over ',num2str(rpt),' runs']);xlabel('Dimension');
% grid on;title('bit error');plotS(errorbits',{'fast','fast+SD'},rk);
end
