
% illustrate how to design the roof of Gaomi Railway Station
% chenzhe, 2017-12-07

W = 140;
T = 51.9;

w = W/2;
t = T/2;

ru = 2450.5;
rd = 1226.3;

R = sqrt(t^2 + ru^2)
X = sqrt(R^2 - w^2)

r = sqrt(t^2 + rd^2)
x = sqrt(r^2 - w^2)

sqrt(ru^2-w^2) + (23.1-21.9) - sqrt(rd^2-w^2)

sqrt(X^2-t^2) + (23.1-21.9) - sqrt(x^2-t^2)

%%
[x,y,z] = sphere(500);
[xx,yy,zz] = sphere(500);

x = x*R;
y = y*R;
z = z*R;

xx = xx*r;
yy = yy*r;
zz = zz*r;
figure;
surf(x,y,z,'facealpha',0.2)
hold on; axis equal;
surf(xx,yy,zz-1226.39,'facealpha',0.8)
view(0,0);set(gca,'ylim',[20,30])

[a,b,c] = meshgrid([-w:10:w],t,[0:-10:-3000]);
a = squeeze(a);
b = squeeze(b);
c = squeeze(c);
surf(a,b,c)