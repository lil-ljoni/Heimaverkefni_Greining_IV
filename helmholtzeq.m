function c=helmholtzeq(L1,L2,h,lambda);

xb=0;
yb=0;

N=L1/h;
M=L2/h;
q=-lambda^2;


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
    A(vinstri_p,vinstri_p)=4/h^2+q;
    %A(vinstri_p,vinstri_p-1)=-2/h^2;%enginn punktur til vinstri
    A(vinstri_p,vinstri_p+1)=-2/h^2;
    A(vinstri_p,vinstri_p-N-1)=-1/h^2;
    A(vinstri_p,vinstri_p+N+1)=-1/h^2;
    %y=yb+(k-1)*h;
    b_vigur(vinstri_p)=ff(0,yb+(k-1)*h);%gamma og beta eru = 0
end

%efri jadar
for j=1:N+1
    efri_p=j+M*(N+1);
    x=xb+(j-1)*h;
    A(efri_p,efri_p)=1;
    b_vigur(efri_p)=vf(xb+h*(j-1));
end

%haegri jadar
for k=2:M
    %haegri_p=M+1+(k-1)*(M+1);
    haegri_p=k*(N+1);
    A(haegri_p,haegri_p)=4/h^2+q;
    A(haegri_p,haegri_p-1)=-2/h^2;
    %A(haegri_p,haegri_p+1)=0;%engin punktur til vinstri
    A(haegri_p,haegri_p-N-1)=-1/h^2;
    A(haegri_p,haegri_p+N+1)=-1/h^2;
    %y=yb+(k-1)*h;
    b_vigur(haegri_p)=ff(L1,yb+(k-1)*h);%gamma og beta eru = 0
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

A = sparse(A);
c=A\b_vigur;

%HZ=zeros(M+1,N+1);

%for j=1:N+1
%    for k=1:M+1
%        HZ(j,k)=c(k+(j-1)*(M+1));
%    end
%end

%HZ;
x=(0:N+1)*h;
y=(0:M+1)*h;
mn=(M+1)*(N+1);
HZ=(reshape(c(1:mn),N+1,M+1))'
%mesh(x,y,HZ)
