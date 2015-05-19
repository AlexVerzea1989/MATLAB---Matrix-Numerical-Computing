function v=sphere_volume(n)
n=max(n,2);
v=zeros(n,1);
v(1)=2;
v(2)=pi;
for i=3:n
    v(i)=2*pi/i*v(i-2);
end
