function thesis_test_eta(rho)

load thesis_example1
k=length(rho);
results=zeros(3,k);
for i=1:k
    results(1,i)=scmplx_beta2(R1,rho(i));
    results(2,i)=scmplx_beta2(R2,rho(i));
    results(3,i)=scmplx_beta2(hR2,rho(i));
end

[rho;results]
