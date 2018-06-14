
global mark1 iMark1
mark1 = mark1;
iMark1 = 14;

tx = -88.58*1e-3;
ty = 0;
tz =  27.98*1e-3;
qw = 0;
qx = 0;
qy = 0;
qz = 0;
x0 = [tx, ty, tz, qw, qx, qy, qz ];


%[x val]=gradient_descent(@(x)func(x),x0,0.2, 1000, 1e-10);
x = x0;
gamma = 0.2;
disp = 0;

for i=1:10000
    i
    [val, grad]=func(x);
    val
    xnew=val-gamma*grad;		%descend along gradient
    change=norm(xnew-x,2);	%calculate how much the x value has changed
    x=xnew;
    T = getTranslationMetrix(x(1), x(2), x(3));
    R = eye(4); Rtmp = qGetR(x(4:7)); R(1:3,1:3) = Rtmp;
    Hc = T*R
    %if display
        % fprintf('%i\t%.3f\t%.3f\t%.3f\n',i,val,change,norm(grad,2));
    %end
    if change<1e-8		%if change is really small
        fprintf('Stopped since ||x-xold|| is less than threshold\n');
        break;			%stop
    end
end


T = getTranslationMetrix(x(1), x(2), x(3));
R = eye(4); Rtmp = qGetR(x(4:7)); R(1:3,1:3) = Rtmp;
Hc = T*R;

subplot(1,2,1);
for i = 1:14
    p0tmp = mark1.pos(:,i);
    p0(1) = -p0tmp(3);
    p0(2) =  p0tmp(1);
    p0(3) =  p0tmp(2);
    p0(4) = 1;   
    H  = mark1.frames(:,:,i);
    p = H*p0';
    plot3(p(1),p(2),p(3),'*');
    hold on
end
grid on
hold off

subplot(1,2,2);
for i = 1:14
    p0tmp = mark1.pos(:,i);
    p0(1,1) = -p0tmp(3);
    p0(1,2) =  p0tmp(1);
    p0(1,3) =  p0tmp(2);
    p0(1,4) = 1;
    H  = mark1.frames(:,:,i);
    p = H*Hc*p0';
    plot3(p(1),p(2),p(3),'*');
    hold on
end
grid on
hold off

%%_________________________________________________________________________
function [val, grad] = func(x)
val = f(x);
grad = df(x);
end

%%_________________________________________________________________________
function val = f(x)
global mark1 iMark1

[T, ~] = getTranslationMetrix(x(1), x(2), x(3));
R = eye(4); Rtmp = qGetR(x(4:7)); R(1:3,1:3) = Rtmp;
H = T*R;
p = zeros(4, iMark1);

dist = zeros(1,iMark1);
for i=1:iMark1
    p (:,i) = mark1.frames(:,:,i)*H*mark1.pos(:,i);
    dist(i) = sqrt(p(1,i).^2+p(2,i).^2+p(3,i).^2);
end
val = std(dist);
end


%%_________________________________________________________________________
function grad = df(x)
global mark1 iMark1

dx = 1e-6;

T = getTranslationMetrix(x(1), x(2), x(3));
R = eye(4); Rtmp = qGetR(x(4:7)); R(1:3,1:3) = Rtmp;
H = T*R;
p = zeros(4, iMark1);
dist = zeros(1,iMark1);
for i=1:iMark1
    p (:,i) = mark1.frames(:,:,i)*H*mark1.pos(:,i);
    dist(i) = sqrt(p(1,i).^2+p(2,i).^2+p(3,i).^2);
end
y = std(dist);

T = getTranslationMetrix(x(1)+dx, x(2), x(3));
R = eye(4); Rtmp = qGetR(x(4:7)); R(1:3,1:3) = Rtmp;
H = T*R;
p = zeros(4, iMark1);
dist = zeros(1,iMark1);
for i=1:iMark1
    p (:,i) = mark1.frames(:,:,i)*H*mark1.pos(:,i);
    dist(i) = sqrt(p(1,i).^2+p(2,i).^2+p(3,i).^2);
end
y1 = std(dist); grad(1) = (y1-y)/dx;

T = getTranslationMetrix(x(1), x(2)+dx, x(3));
R = eye(4); Rtmp = qGetR(x(4:7)); R(1:3,1:3) = Rtmp;
H = T*R;
p = zeros(4, iMark1);
dist = zeros(1,iMark1);
for i=1:iMark1
    p (:,i) = mark1.frames(:,:,i)*H*mark1.pos(:,i);
    dist(i) = sqrt(p(1,i).^2+p(2,i).^2+p(3,i).^2);
end
y2 = std(dist); grad(2) = (y2-y)/dx;

T = getTranslationMetrix(x(1), x(2), x(3)+dx);
R = eye(4); Rtmp = qGetR(x(4:7)); R(1:3,1:3) = Rtmp;
H = T*R;
p = zeros(4, iMark1);
dist = zeros(1,iMark1);
for i=1:iMark1
    p (:,i) = mark1.frames(:,:,i)*H*mark1.pos(:,i);
    dist(i) = sqrt(p(1,i).^2+p(2,i).^2+p(3,i).^2);
end
y3 = std(dist); grad(3) = (y3-y)/dx;

T = getTranslationMetrix(x(1), x(2), x(3));
R = eye(4); Rtmp = qGetR([x(4)+dx, x(5), x(6), x(7)]); R(1:3,1:3) = Rtmp;
H = T*R;
p = zeros(4, iMark1);
dist = zeros(1,iMark1);
for i=1:iMark1
    p (:,i) = mark1.frames(:,:,i)*H*mark1.pos(:,i);
    dist(i) = sqrt(p(1,i).^2+p(2,i).^2+p(3,i).^2);
end
y4 = std(dist); grad(4) = (y4-y)/dx;

T = getTranslationMetrix(x(1), x(2), x(3));
R = eye(4); Rtmp = qGetR([x(4), x(5)+dx, x(6), x(7)]); R(1:3,1:3) = Rtmp;
H = T*R;
p = zeros(4, iMark1);
dist = zeros(1,iMark1);
for i=1:iMark1
    p (:,i) = mark1.frames(:,:,i)*H*mark1.pos(:,i);
    dist(i) = sqrt(p(1,i).^2+p(2,i).^2+p(3,i).^2);
end
y5 = std(dist); grad(5) = (y5-y)/dx;

T = getTranslationMetrix(x(1), x(2), x(3));
R = eye(4); Rtmp = qGetR([x(4), x(5), x(6)+dx, x(7)]); R(1:3,1:3) = Rtmp;
H = T*R;
p = zeros(4, iMark1);
dist = zeros(1,iMark1);
for i=1:iMark1
    p (:,i) = mark1.frames(:,:,i)*H*mark1.pos(:,i);
    dist(i) = sqrt(p(1,i).^2+p(2,i).^2+p(3,i).^2);
end
y6 = std(dist); grad(6) = (y6-y)/dx;

T = getTranslationMetrix(x(1), x(2), x(3));
R = eye(4); Rtmp = qGetR([x(4), x(5), x(6), x(7)+dx]); R(1:3,1:3) = Rtmp;
H = T*R;
p = zeros(4, iMark1);
dist = zeros(1,iMark1);
for i=1:iMark1
    p (:,i) = mark1.frames(:,:,i)*H*mark1.pos(:,i);
    dist(i) = sqrt(p(1,i).^2+p(2,i).^2+p(3,i).^2);
end
y7 = std(dist); grad(7) = (y7-y)/dx;

end


