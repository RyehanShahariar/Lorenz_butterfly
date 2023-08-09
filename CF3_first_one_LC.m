clc;
clear;
close all;
%inputs
h=0.01; t(1)=0; tfinal=500;
t=t(1):h:tfinal; N=ceil((tfinal-t(1))/h);
%inital
x(1)=0.1; y(1)=0.1; z(1)=0.1;
%constant
%constant
alpha=1;
%parameters
sigma=10; b=8/3; r=30;
%Given ODE 3.1-3.1c
f1=@(t,x,y,z) -sigma*(x-y);
f2=@(t,x,y,z) r.*x-y-x.*z;
f3=@(t,x,y,z) x.*y-b.*z;
% %constant Atangana-Baleanu-Caputo(mittag-Leffler)

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
figure(1);
plot3(x,y,z,'m');
% c = 1:numel(t);      %# colors
% h = surface([x(:), x(:)], [y(:), y(:)], [z(:), z(:)], ...
%     [c(:), c(:)], 'EdgeColor','flat', 'FaceColor','none');
% colormap( jet(numel(t)) )
xlabel('x'),ylabel('y'),zlabel('z')
% legend('contant fractional-order')
toc;

% 
% %variable version  caputo-Fabrizio 
% alpha = @(t) 1./(1+exp(-t));
% 
% for n=2:N
%     x(n+1) = x(n)+(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*3*h*alpha(t(n))*(2-alpha(t(n))))*...
%         f1(t(n),x(n),y(n),z(n))-(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*h*alpha(t(n))*...
%         (2-alpha(t(n))))*f1(t(n-1),x(n-1),y(n-1),z(n-1));
% 
%     y(n+1) = y(n)+(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*3*h*alpha(t(n))*(2-alpha(t(n))))*...
%         f2(t(n),x(n),y(n),z(n))-(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*h*alpha(t(n))*...
%         (2-alpha(t(n))))*f2(t(n-1),x(n-1),y(n-1),z(n-1));
% 
%     z(n+1) = z(n)+(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*3*h*alpha(t(n))*(2-alpha(t(n))))*...
%         f3(t(n),x(n),y(n),z(n))-(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*h*alpha(t(n))*...
%         (2-alpha(t(n))))*f3(t(n-1),x(n-1),y(n-1),z(n-1));
% 
% t(n+1)=t(n)+h;
% end
% figure(2);
% plot3(x,y,z);
% xlabel('x'),ylabel('y'),zlabel('z'),legend('contant fractional-order')