# Lorenz_butterfly
  Its a lorenz system . runed with matlab. https://doi.org/10.48550/arXiv.2307.03251
  If you run this MATLAB Code you will get the Lorenz Butterfly effect
  Explanation:


h=0.01; t(1)=0; tfinal=500;
t=t(1):h:tfinal; N=ceil((tfinal-t(1))/h);
x(1)=0.1; y(1)=0.1; z(1)=0.1;
alpha=1;sigma=10; b=8/3; r=30;


These 3 are the main equations with 3 variables.

f1=@(t,x,y,z) -sigma*(x-y);
f2=@(t,x,y,z) r.*x-y-x.*z;
f3=@(t,x,y,z) x.*y-b.*z;



x(2)=x(1)+((h^alpha)/(gamma(alpha+1))).*f1(t(1),x(1),y(1),z(1));

y(2)=y(1)+((h^alpha)/(gamma(alpha+1))).*f2(t(1),x(1),y(1),z(1));

z(2)=z(1)+((h^alpha)/(gamma(alpha+1))).*f3(t(1),x(1),y(1),z(1));


for n=2:N
    tic;

x(n+1)=x(n)+(0.5*(2-alpha)*(1-alpha)+0.25*3*h*alpha*(2-alpha))*f1(t(n),x(n),y(n),z(n))-...
    (0.5*(2-alpha)*(1-alpha)+0.25*h*alpha*(2-alpha))*f1(t(n-1),x(n-1),y(n-1),z(n-1));

y(n+1)=y(n)+(0.5*(2-alpha)*(1-alpha)+0.25*3*h*alpha*(2-alpha))*f2(t(n),x(n),y(n),z(n))-...
    (0.5*(2-alpha)*(1-alpha)+0.25*h*alpha*(2-alpha))*f2(t(n-1),x(n-1),y(n-1),z(n-1));

z(n+1)=z(n)+(0.5*(2-alpha)*(1-alpha)+0.25*3*h*alpha*(2-alpha))*f3(t(n),x(n),y(n),z(n))-...
    (0.5*(2-alpha)*(1-alpha)+0.25*h*alpha*(2-alpha))*f3(t(n-1),x(n-1),y(n-1),z(n-1));

t(n+1)=t(n)+h;
end

Now if you like to plot z then you have to make it as the following:
figure(1);
plot3(x,y,z,'m');

xlabel('x'),ylabel('y'),zlabel('z')

toc;
