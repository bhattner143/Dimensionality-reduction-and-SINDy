function ut =FiniteDiff2(u,h,d)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
m=size(u,1);
ut=zeros(m,1);

if d==1
    for i=2:m-1
        ut(i)=(u(i+1)-u(i-1))/(2*h);
    end
    ut(1)=((-3/2)*u(1)+2*u(2)-u(3)/2)/h;
    ut(m)=(1.5*u(m-1)+2*u(m-2)+u(m-3)/2)/h;
end
end

