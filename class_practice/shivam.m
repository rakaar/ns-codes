% Inputs
h=0.1;b=1;a=0; yn=1;

% Generate A B 
n=(b-a)/h;
x=[a:h:b]; 
A=zeros(n,n);
for i=1:n
    A(i,i)=(-2/(h*h))-x(i);
    if(i~=1)
        A(i,i-1)=(1/(h*h));
    end
    if(i~=n)
        A(i,i+1)=(1/(h*h));
    end
end
A(1,1)=(-2/(h*h))-x(1)+2/h;
A(1,2)=(2/(h*h));
B=zeros(n,1);
for i=1:n
    B(i,1)=1;
    if(i==1)
        B(i,1)=1+2/h;
    end
    if(i==n)
        B(i,1)=1-(yn/(h*h));
    end
end

% Result and plot
Y=gauss_elimination(A,B,n);
result=zeros(n+1,1);
result(n+1,1)=yn;
for i=1:n
     result(i,1)=Y(i);
end
disp(result);
figure(1)
    plot(x,result); 
grid

% Gauss Elimination
function x=gauss_elimination(A,B,n)
    %Decomposition
    Ab=[A,B];
    for i=1:size(Ab,1)
        for j=i+1:size(Ab,1)
            key1=Ab(j,i)./Ab(i,i);
            Ab(j,:)=Ab(j,:)-key1.*Ab(i,:);
        end
    end
    x=zeros(1,size(A,2));
    for i=size(Ab,1):-1:1
        hg=sum(Ab(i,i+1:end-1).*x(i+1:end));
        x(i)=(Ab(i,end)-hg)./Ab(i,i);
    end
end