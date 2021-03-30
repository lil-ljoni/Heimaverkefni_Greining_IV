function c=helmholtzeq(L1,L2,h,lambda)

xb=0;
yb=0;

N=L1/h;
M=L2/h;
q=-lambda^2;
wf = @(t) 1;
vf = @(t) 0;
ff = @(t,r) 0;

u0= 10; u1=1;

%wf = @(t) -u0*t/L1*(t/L1-1)^2*(1+t./L1);
%vf = @(t) u1.*t./L1*(1-t./L1)*(1+t./L1).^2;


%x=zeros(N+1,1);
%y=zeros(M+1,1);
b_vigur=zeros((N+1)*(M+1),1);
%b_vigur=zeros(M+1,1);
A=zeros((N+1)*(M+1),(N+1)*(M+1));
%A=zeros(M+1,N+1)

%nedri jadar
for j=1:N+1
    nedri_p=j;
    A(nedri_p,nedri_p)=1;
    x=xb+(j-1)*h;
    b_vigur(nedri_p)=wf(x);
end

%vinstri jadar
for k=2:M
    vinstri_p=1+(k-1)*(N+1);
    A(vinstri_p,vinstri_p)=4/h^2+q;  %j = punktur sem við erum í
    %A(vinstri_p,vinstri_p-1)=-2/h^2;%enginn punktur til vinstri
    A(vinstri_p,vinstri_p+1)=-2/h^2;        %r = hægra megin
    A(vinstri_p,vinstri_p-N-1)=-1/h^2;      %s - neðri punktur
    A(vinstri_p,vinstri_p+N+1)=-1/h^2;      % t - efri punktur
    %y=yb+(k-1)*h;
    b_vigur(vinstri_p)=ff(0,yb+(k-1)*h);%gamma og alpha eru = 0
end

%efri jadar
for j=1:N+1
    efri_p=j+M*(N+1);
    x=xb+(j-1)*h;
    A(efri_p,efri_p)=1;
    b_vigur(efri_p)=vf(x);
end

%haegri jadar
for k=2:M
    %haegri_p=M+1+(k-1)*(M+1);
    haegri_p=k*(N+1);
    A(haegri_p,haegri_p)=4/h^2+q;
    A(haegri_p,haegri_p-1)=-2/h^2;
    %A(haegri_p,haegri_p+1)=0;%engin punktur til hægri
    A(haegri_p,haegri_p-N-1)=-1/h^2;    %s
    A(haegri_p,haegri_p+N+1)=-1/h^2;    %t
    %y=yb+(k-1)*h;
    b_vigur(haegri_p)=ff(L1,yb+(k-1)*h);%gamma og alpha eru = 0
end

%innra svaedi
for k=2:M
    for j=2:N
        innri_p=j+(k-1)*(N+1);
        %x=xb+(j-1)*h;
        %y=yb+(k-1)*h;
        A(innri_p,innri_p)=4/h^2+q;
        A(innri_p,innri_p-1)=-1/h^2;
        A(innri_p,innri_p+1)=-1/h^2;
        A(innri_p,innri_p-N-1)=-1/h^2;
        A(innri_p,innri_p+N+1)=-1/h^2;
        b_vigur(innri_p)=ff(xb+(j-1)*h,yb+(k-1)*h);
    end
end
format short
A = sparse(A);
c=A\b_vigur;



mn=(M+1)*(N+1);
HZ=(reshape(c(1:mn),N+1,M+1))' % nálgun á lausn

%velja plott

plotta_1(HZ,L1,L2,lambda,N,M,h)
%plotta_2(HZ,L1,L2,lambda,N,M,h)

function plotta_1(HZ,L1,L2,lambda,N,M,h)

[x,y] = meshgrid((0:N).*h,(0:M).*h);
u_l = sin(lambda.*(L2-y))./sin(L2*lambda);
figure(1)
surf(x,y,u_l,"Facecolor","g")
hold on
surf(x,y,HZ,"FaceColor","r")
hold off
grid on; xlabel("x");ylabel("y");zlabel("z")
title("Lausnin {u_e} og  nálgunin fyrir {\lambda}="+lambda)
legend("{u_e}","{c_{jk}}")
end
end

function plotta_2(HZ,L1,L2,lambda,N,M,h)
[x,y] = meshgrid((0:N).*h,(0:M).*h);
figure;
surf(x,y,HZ)
grid on; xlabel("x");ylabel("y");zlabel("z")
title("Nálgun á lausn fyrir {\lambda}="+lambda)
colormap hsv
end
