function testSearch(rk,rpt,econd,noise,id)
tests={@testSDShrink,@testSDFixed,@scmplx};
cmplx=zeros(length(rk),length(tests));
nt=numel(tests);
for n=1:length(rk)
    rn=rk(n);
    disp(rn)
    for j=1:rpt
        [A,y,x]=randILS(rn,rn,econd,noise);
        Ay=PLLL([A,y],rn,1,@QRZ_PHH,@permu_Givens,2);
        parfor i=1:nt
            testi=tests{i};
            %-------complexity----
            [~,~,c]=testi(Ay);
            cmplx(n,i)=cmplx(n,i)+c;
        end
    end
    cmplx(n,:)=cmplx(n,:)/rpt;
end

[xaxis,ix]=sort(rk);
names=cell(size(tests));
for i=1:numel(tests)
    names{i}=func2str(tests{i});
    if strcmp(names{i}(1:4),'test')
        names{i}=names{i}(5:length(names{i}));
    end
end

figure;hold on;ylabel('Average complexity');xlabel('Dimension');grid on;plotS(cmplx(ix,:),names,xaxis);

end

function [Z,rsd,cmplx]=testSDShrink(Ry)
    init_count(0);
    init_count2(0);
    [Z,rsd]=search_CSD(Ry,inf);
    init_count(0);
    cmplx=init_count2(0);
end

function [Z,rsd,cmplx]=testSDFixed(Ry)
    init_count(0);
    init_count2(0);
    [Z,rsd]=search_CSD(Ry,-inf);
    init_count(0);
    cmplx=init_count2(0);
end
