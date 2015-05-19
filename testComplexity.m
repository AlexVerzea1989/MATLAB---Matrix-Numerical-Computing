function E=testComplexity(rk,noises,rpt,econd,show_dimension,id)
if nargin<5
    show_dimension=true;
end
E=zeros(2,length(rk));
for n=1:length(rk)
    rn=rk(n);
    noise=noises(n);
    disp([rn,noise])
    for j=1:rpt
%         disp(j)
        [A,y,x]=randILS(rn,rn,econd,noise);
        Ay=[A,y];
        Ry1=QRZ_HH(Ay,rn);
        
        [Ry2,Z2]=PLLL(Ay,rn);
        [z,beta,E1]=scmplx(Ry2);
        init_count(0);
        init_count2(0);
        Zhat = search_CSD(Ry2,-beta);
        E2=init_count2();
        E(1,n)=E(1,n)+E1;
        E(2,n)=E(2,n)+E2;
    end
end
E=E/rpt;
figure
if show_dimension
    plot(rk,E(1,:),'r-o',rk,E(2,:),'b:.');
    xlabel(['dimension (sigma = ',num2str(noises(1)),')']);
else
    plot(noises,E(1,:),'r-o',noises,E(2,:),'b.:');
    xlabel(['Sigma (Dimension = ',num2str(rk(1)),')']);
end
ylabel('Search Complexity');
set(gca,'yscale','log');
legend({'theoritical','experimental'},'Location','Best');
