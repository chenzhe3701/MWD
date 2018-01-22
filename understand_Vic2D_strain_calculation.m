% this code tries to investigate how Vic2D calculates Lagrangian strain
%
% chenzhe, 2017-10-02
%
% chenzhe, 2018-01-10
% maybe there was an error?  plot point position and edit.

clc; close all;
% (1) Pick 3 points, and show how strain was calculated using these 3 points
% edge length of original triangle
dL = x(1,2) - x(1,1)
u(sigma==-1) = nan;
v(sigma==-1) = nan;

% select a point
iR = 10;
iC = 10;

% plot points of interest.
figure;
hold on;
plot([x(iR-1,iC),x(iR,iC),x(iR,iC+1)],[y(iR-1,iC),y(iR,iC),y(iR,iC+1)],'-or');
plot([x(iR-1,iC),x(iR,iC),x(iR,iC+1)]+[u(iR-1,iC),u(iR,iC),u(iR,iC+1)] - u(iR,iC),...
    [y(iR-1,iC),y(iR,iC),y(iR,iC+1)]+[v(iR-1,iC),v(iR,iC),v(iR,iC+1)] - v(iR,iC),'-ob');
set(gca,'ydir','reverse');

% right point
du1 = u(iR,iC+1) - u(iR,iC)
dv1 = v(iR,iC+1) - v(iR,iC)
% up point
du2 = u(iR-1,iC) - u(iR,iC)
dv2 = v(iR-1,iC) - v(iR,iC)

% use right point to calculate d?/dX
dudX = du1/dL
dvdX = dv1/dL
% use up point to calcualte d?/dY
dudY = du2/(-dL)    % when going up, Y decreases, so dY<0, dY = (-dL)
dvdY = dv2/(-dL)

disp('displacement gradient Fd = ')
Fd = [dudX dudY; dvdX dvdY]     % 'displacement gradient'
disp('note: Fd, Fu and Fu_d are all displacement gradient rather than deformation gradient.')
disp('Only deformation gradient should use the letter F')
%% (1a) old (2017-10 method of rotation back).
disp('2017-10 method: on 2018-10 I found this was likely a mistake, but still did not affect the finite strain --------------------------------')
theta = atand(dv1/du1)
M = [cosd(theta), cosd(theta-90);
    cosd(theta+90), cosd(theta)]
R = M'

t = M*[du1;dv1]
du3 = t(1)
dv3 = t(2)

t = M*[du2;dv2]
du4 = t(1)
dv4 = t(2)
    
dudX = du3/dL
dvdX = dv3/dL
dudY = du4/(-dL)
dvdY = dv4/(-dL)

% something like the 'stretch tensor', but
% (1) this Fu is not symmetric, while 'stretch tensor' seems have to be symmetric according to the definition I know.
% (2) this Fu is also not the 'raw'.  The rotation is taken off, but in an 'artificial' way, i.e., forcing deformed x-edge to be horizontal.  
Fu = [dudX dudY; dvdX dvdY]  
disp('R*Fu-Fd: ')
R*Fu - Fd   % This is equal.  I don't know why eliminate rotation first.
% chenzhe, 2018-01: R*Fu ~= Fd.  But there was a mistake before.  The rotation angle I used was not correct.
% Moreover, using F=I+Fu vs F=I+Fd, the finite strain was not equal.

% finite strain assumption
F = Fd+eye(2)
E_finite_raw = (F'*F-eye(2))/2

F = Fu+eye(2)
E_finite_rotate_first = (F'*F-eye(2))/2

disp('E_finite_raw - E_finite_rotate_first = ')
E_finite_raw - E_finite_rotate_first

% small strain assumption
disp('small strain, raw:')
(Fd'+Fd)/2
disp('small strain, rotate first:')
(Fu'+Fu)/2


%% (1b) new (2018-01) method of rotation back
disp('2018-01 method: --------------------------------')
theta = atand(dv1/(du1+dL))
M = [cosd(theta), cosd(theta-90);
    cosd(theta+90), cosd(theta)]
R = M'

t = M*[dL+du1;dv1]
du5 = t(1)-dL
dv5 = t(2)

t = M*[du2;-dL+dv2]
du6 = t(1)
dv6 = t(2) - (-dL)
    
dudX = du5/dL
dvdX = dv5/dL
dudY = du6/(-dL)
dvdY = dv6/(-dL)

Fu_b = [dudX dudY; dvdX dvdY]
disp('R*Fu_b - Fd: ')
R*Fu_b - Fd   % chenzhe, 2018-01.  R*Fu_b != Fd.  However, using finite strain definition, the strain calculated by Fu_b = the strain calculated by Fd.

% finite strain assumption
F = Fd+eye(2)
E_finite_raw = (F'*F-eye(2))/2

F = Fu_b+eye(2)
E_finite_rotate_first = (F'*F-eye(2))/2

disp('E_finite_raw - E_finite_rotate_first = ')
E_finite_raw - E_finite_rotate_first

% small strain assumption
disp('small strain, raw:')
(Fd'+Fd)/2
disp('small strain, rotate first:')
(Fu_b'+Fu_b)/2


%% (2) Find out how Vic2D's filter looks ( step 3 was actually done first, but cound't match Vic's description, so here just find out what Vic is using)
% exx1 is some raw data.  exx2 is the data after applying Vic's 'decay filter'. 
% "... a decay filter, worth 10% at the edges" 
% "The decay filter is a 90% center-weighted Gaussian filter and works best for most situations" 
% I tried at least two fitler_size (5 and 11)

load('WE43_T6_C1_s5_r0c0_exx_fs_5.mat','exx1','exx2','sigma','filter_size');
hfs = (filter_size-1)/2;
s = ones(size(sigma));
s(sigma==-1)=0;
s = filter2(fspecial('average',filter_size),s,'same');
s(s<0.999999) = -1;
n = sum(sum(s>-1));
% y = x*b
y=zeros(n,1);
x=zeros(n,filter_size*filter_size);
[nR,nC] = size(sigma);
ii = 1;
for iR = 1:nR
   for iC = 1:nC
      if s(iR,iC)>-1
         yLocal = exx2(iR,iC);
         xLocal = exx1(iR-hfs:iR+hfs,iC-hfs:iC+hfs);
         y(ii) = yLocal;
         x(ii,:) = xLocal(:)';
         ii = ii+1;
      end
   end
   disp(iR);
end
b = x\y
b = reshape(b,filter_size,filter_size)      % this is the actual filter that Vic2D used.


filter_size = filter_size;
hfs = (filter_size-1)/2;     % half filter size
tos  = icdf('normal',0.5 + 0.9/2, 0, 1)  % ~1.645(fs=11, fs=5) Times of sigma corresponding to cummulative prob = X, so that 2(X-0.5)=0.9.  This is 1-d, does not seem good
tos  = icdf('normal',0.5 + sqrt(0.9/4), 0, 1)  % ~1.949(fs=11, fs=5)  Times of sigma corresponding to cummulative prob = X, so that 4*(x-0.5)^2=0.9  This is 2-d square, also does not work


% solve for tos actually used in Vic2D's filter, numerical solution
xi = 0;
xf = 10;
x = linspace(xi,xf,101);
dy = 1;
target = 1e-6;
while(dy>target)
    y = abs((pdf('normal',x,0,1)./pdf('normal',0,0,1)-b(hfs+1,1)/b(hfs+1,hfs+1)))
    [dy,ind] = min(y)
    xi = x(ind)
    x(ind) = [];
    y(ind) = [];
    [dy2,ind] = min(y)
    xf = x(ind)
    x = linspace(xi,xf,101);
end
tos = xi    % ~1.6113(fs=11) ~1.4179(fs=5)
(cdf('normal',tos,0,1) - 0.5)*2
pdf('normal',0,0,1)/pdf('normal',tos,0,1)

s_dp = hfs/tos     % when hFilterSize corresponds to tos, sigma value in unit of # of data points

h = fspecial('gaussian',filter_size,s_dp)    % the filter with the calculated sigma, and with 90% cummulative prob within area of filterSize.  Filter is square
filter_size/s_dp        % numerically, filter_size = 3.5448 times of s_dp 

%% (3) use matrix method to calculate for all data points (finite strain)
% looks like on 2017-10, I did not do the rotation, based on the fact that without rotation, the finite strain calculation was almost the same.  

load('WE43_T6_C1_s5_r0c0.mat')
dL = x(1,2) - x(1,1);

du1 = circshift(u,[0,-1])-u;
dv1 = circshift(v,[0,-1])-v;
du2 = circshift(u,[1,0])-u;
dv2 = circshift(v,[1,0])-v;
du3 = circshift(u,[0,1])-u;
dv3 = circshift(v,[0,1])-v;
du4 = circshift(u,[-1,0])-u;
dv4 = circshift(v,[-1,0])-v;

dudX1 = du1/dL; dudX1(:,end) = nan;
dvdX1 = dv1/dL; dvdX1(:,end) = nan;
dudY2 = du2/(-dL); dudY2(1,:) = nan;
dvdY2 = dv2/(-dL); dvdY2(1,:) = nan;
dudX3 = du3/(-dL); dudX3(:,1) = nan;
dvdX3 = dv3/(-dL); dvdX3(:,1) = nan;
dudY4 = du4/dL; dudY4(end,:) = nan;
dvdY4 = dv4/dL; dvdY4(end,:) = nan;

F1 = arrayfun(@(a,b,c,d) [a,b;c,d]+eye(2), dudX1,dudY2,dvdX1,dvdY2,'uniformoutput',0);
e1 = cellfun(@(x) (x'*x-eye(2))/2, F1,'uniformoutput',0);
exx1 = cell2mat(cellfun(@(x) x(1), e1, 'uniformoutput',0));
exy1 = cell2mat(cellfun(@(x) x(2), e1, 'uniformoutput',0));
eyy1 = cell2mat(cellfun(@(x) x(4), e1, 'uniformoutput',0));


F2 = arrayfun(@(a,b,c,d) [a,b;c,d]+eye(2), dudX3,dudY2,dvdX3,dvdY2,'uniformoutput',0);
e2 = cellfun(@(x) (x'*x-eye(2))/2, F2,'uniformoutput',0);
exx2 = cell2mat(cellfun(@(x) x(1), e2, 'uniformoutput',0));
exy2 = cell2mat(cellfun(@(x) x(2), e2, 'uniformoutput',0));
eyy2 = cell2mat(cellfun(@(x) x(4), e2, 'uniformoutput',0));

F3 = arrayfun(@(a,b,c,d) [a,b;c,d]+eye(2), dudX3,dudY4,dvdX3,dvdY4,'uniformoutput',0);
e3 = cellfun(@(x) (x'*x-eye(2))/2, F3,'uniformoutput',0);
exx3 = cell2mat(cellfun(@(x) x(1), e3, 'uniformoutput',0));
exy3 = cell2mat(cellfun(@(x) x(2), e3, 'uniformoutput',0));
eyy3 = cell2mat(cellfun(@(x) x(4), e3, 'uniformoutput',0));

F4 = arrayfun(@(a,b,c,d) [a,b;c,d]+eye(2), dudX1,dudY4,dvdX1,dvdY4,'uniformoutput',0);
e4 = cellfun(@(x) (x'*x-eye(2))/2, F4,'uniformoutput',0);
exx4 = cell2mat(cellfun(@(x) x(1), e4, 'uniformoutput',0));
exy4 = cell2mat(cellfun(@(x) x(2), e4, 'uniformoutput',0));
eyy4 = cell2mat(cellfun(@(x) x(4), e4, 'uniformoutput',0));

exx0 = mean(cat(3,exx1,exx2,exx3,exx4),3);
exy0 = mean(cat(3,exy1,exy2,exy3,exy4),3);
eyy0 = mean(cat(3,eyy1,eyy2,eyy3,eyy4),3);

filter_size = 5;
hfs = (filter_size-1)/2;     % half filter size
tos  = icdf('normal',0.5 + 0.9/2, 0, 1)  % ~1.645(fs=11, fs=5) Times of sigma corresponding to cummulative prob = X, so that 2(X-0.5)=0.9.  This is 1-d, does not seem good
tos  = icdf('normal',0.5 + sqrt(0.9/4), 0, 1)  % ~1.949(fs=11, fs=5)  Times of sigma corresponding to cummulative prob = X, so that 4*(x-0.5)^2=0.9  This is 2-d square, also does not work

s_dp = hfs/tos     % when hFilterSize correspond to tos, sigma value in unit of # of data points
h = fspecial('gaussian',filter_size,s_dp)    % the filter with the calculated sigma, and with 90% cummulative prob within area of filterSize.  Filter is square

hd = logical(fspecial('disk',hfs));
h2 = h.*hd;
h2 = h2/sum(h2(:))                          % Similar to h, but for h2, filter size is a disk.

h = [
    0.0124    0.0263    0.0339    0.0263    0.0124
    0.0263    0.0560    0.0720    0.0560    0.0263
    0.0339    0.0720    0.0925    0.0720    0.0339
    0.0263    0.0560    0.0720    0.0560    0.0263
    0.0124    0.0263    0.0339    0.0263    0.0124]         % this is copied from fitted coefficients
exxA = filter2(h,exx0);

myplot(exx);caxis([-0.12 0.02])     % (1) raw out put from Vic2D
myplot(exx0);caxis([-0.12 0.02])    % (2) 4-point average when reconstructing strain from Vic2D's displacement
myplot(exxA);caxis([-0.12 0.02])    % (3) 4-point averaged, then filtered strain, reconstructed from Vic2D's displacement.
myplot(abs(exx-exx0) < 1e-3)        %
myplot(abs(exx-exxA) < 1e-3)        % Reconstruction from (3) is very similar to direct output (1) from Vic2D



%% (4) use matrix method to calculate for all data points (infinitesimal strain).  Looks like this is not as good as finite strain.
load('WE43_T6_C1_s5_r0c0.mat')
dL = x(1,2) - x(1,1);

du1 = circshift(u,[0,-1])-u;
dv1 = circshift(v,[0,-1])-v;
du2 = circshift(u,[1,0])-u;
dv2 = circshift(v,[1,0])-v;
du3 = circshift(u,[0,1])-u;
dv3 = circshift(v,[0,1])-v;
du4 = circshift(u,[-1,0])-u;
dv4 = circshift(v,[-1,0])-v;

theta = arrayfun(@(a,b) atand(a/b), dv1, du1,'uniformoutput',0);
M = cellfun(@(x) [cosd(x), cosd(x-90);cosd(x+90), cosd(x)], theta, 'uniformoutput',0);
t = cellfun(@(x,y,z) x*[y;z], M,num2cell(du1),num2cell(dv1),'uniformoutput',0);
du_1 = cell2mat(cellfun(@(x) x(1),t,'uniformoutput',0));
dv_1 = cell2mat(cellfun(@(x) x(2),t,'uniformoutput',0));
t = cellfun(@(x,y,z) x*[y;z], M,num2cell(du2),num2cell(dv2),'uniformoutput',0);
du_2 = cell2mat(cellfun(@(x) x(1),t,'uniformoutput',0));
dv_2 = cell2mat(cellfun(@(x) x(2),t,'uniformoutput',0));
dudX1 = du_1/(dL); dudX1(:,end) = nan;
dvdX1 = dv_1/(dL); dvdX1(:,end) = nan;
dudY2 = du_2/(-dL); dudY2(1,:) = nan;
dvdY2 = dv_2/(-dL); dvdY2(1,:) = nan;
F1 = arrayfun(@(a,b,c,d) [a,b;c,d], dudX1,dudY2,dvdX1,dvdY2,'uniformoutput',0);
e1 = cellfun(@(x) (x'+x)/2, F1,'uniformoutput',0);
exx1 = cell2mat(cellfun(@(x) x(1), e1, 'uniformoutput',0));
exy1 = cell2mat(cellfun(@(x) x(2), e1, 'uniformoutput',0));
eyy1 = cell2mat(cellfun(@(x) x(4), e1, 'uniformoutput',0));

theta = arrayfun(@(a,b) atand(a/b), dv3, du3,'uniformoutput',0);
M = cellfun(@(x) [cosd(x), cosd(x-90);cosd(x+90), cosd(x)], theta, 'uniformoutput',0);
t = cellfun(@(x,y,z) x*[y;z], M,num2cell(du3),num2cell(dv3),'uniformoutput',0);
du_3 = cell2mat(cellfun(@(x) x(1),t,'uniformoutput',0));
dv_3 = cell2mat(cellfun(@(x) x(2),t,'uniformoutput',0));
t = cellfun(@(x,y,z) x*[y;z], M,num2cell(du2),num2cell(dv2),'uniformoutput',0);
du_2 = cell2mat(cellfun(@(x) x(1),t,'uniformoutput',0));
dv_2 = cell2mat(cellfun(@(x) x(2),t,'uniformoutput',0));
dudX3 = du_3/(-dL); dudX3(:,1) = nan;
dvdX3 = dv_3/(-dL); dvdX3(:,1) = nan;
dudY2 = du_2/(-dL); dudY2(1,:) = nan;
dvdY2 = dv_2/(-dL); dvdY2(1,:) = nan;
F2 = arrayfun(@(a,b,c,d) [a,b;c,d], dudX3,dudY2,dvdX3,dvdY2,'uniformoutput',0);
e2 = cellfun(@(x) (x'+x)/2, F2,'uniformoutput',0);
exx2 = cell2mat(cellfun(@(x) x(1), e2, 'uniformoutput',0));
exy2 = cell2mat(cellfun(@(x) x(2), e2, 'uniformoutput',0));
eyy2 = cell2mat(cellfun(@(x) x(4), e2, 'uniformoutput',0));

theta = arrayfun(@(a,b) atand(a/b), dv3, du3,'uniformoutput',0);
M = cellfun(@(x) [cosd(x), cosd(x-90);cosd(x+90), cosd(x)], theta, 'uniformoutput',0);
t = cellfun(@(x,y,z) x*[y;z], M,num2cell(du3),num2cell(dv3),'uniformoutput',0);
du_3 = cell2mat(cellfun(@(x) x(1),t,'uniformoutput',0));
dv_3 = cell2mat(cellfun(@(x) x(2),t,'uniformoutput',0));
t = cellfun(@(x,y,z) x*[y;z], M,num2cell(du4),num2cell(dv4),'uniformoutput',0);
du_4 = cell2mat(cellfun(@(x) x(1),t,'uniformoutput',0));
dv_4= cell2mat(cellfun(@(x) x(2),t,'uniformoutput',0));
dudX3 = du_3/(-dL); dudX3(:,1) = nan;
dvdX3 = dv_3/(-dL); dvdX3(:,1) = nan;
dudY4 = du_4/(dL); dudY4(end,:) = nan;
dvdY4 = dv_4/(dL); dvdY4(end,:) = nan;
F3 = arrayfun(@(a,b,c,d) [a,b;c,d], dudX3,dudY4,dvdX3,dvdY4,'uniformoutput',0);
e3 = cellfun(@(x) (x'+x)/2, F3,'uniformoutput',0);
exx3 = cell2mat(cellfun(@(x) x(1), e3, 'uniformoutput',0));
exy3 = cell2mat(cellfun(@(x) x(2), e3, 'uniformoutput',0));
eyy3 = cell2mat(cellfun(@(x) x(4), e3, 'uniformoutput',0));

theta = arrayfun(@(a,b) atand(a/b), dv1, du1,'uniformoutput',0);
M = cellfun(@(x) [cosd(x), cosd(x-90);cosd(x+90), cosd(x)], theta, 'uniformoutput',0);
t = cellfun(@(x,y,z) x*[y;z], M,num2cell(du1),num2cell(dv1),'uniformoutput',0);
du_1 = cell2mat(cellfun(@(x) x(1),t,'uniformoutput',0));
dv_1 = cell2mat(cellfun(@(x) x(2),t,'uniformoutput',0));
t = cellfun(@(x,y,z) x*[y;z], M,num2cell(du4),num2cell(dv4),'uniformoutput',0);
du_4 = cell2mat(cellfun(@(x) x(1),t,'uniformoutput',0));
dv_4 = cell2mat(cellfun(@(x) x(2),t,'uniformoutput',0));
dudX1 = du_1/(dL); dudX1(:,end) = nan;
dvdX1 = dv_1/(dL); dvdX1(:,end) = nan;
dudY4 = du_4/(dL); dudY4(end,:) = nan;
dvdY4 = dv_4/(dL); dvdY4(end,:) = nan;
F4 = arrayfun(@(a,b,c,d) [a,b;c,d], dudX1,dudY4,dvdX1,dvdY4,'uniformoutput',0);
e4 = cellfun(@(x) (x'+x)/2, F4,'uniformoutput',0);
exx4 = cell2mat(cellfun(@(x) x(1), e4, 'uniformoutput',0));
exy4 = cell2mat(cellfun(@(x) x(2), e4, 'uniformoutput',0));
eyy4 = cell2mat(cellfun(@(x) x(4), e4, 'uniformoutput',0));

exx0_small = mean(cat(3,exx1,exx2,exx3,exx4),3);
exy0_small = mean(cat(3,exy1,exy2,exy3,exy4),3);
eyy0_small = mean(cat(3,eyy1,eyy2,eyy3,eyy4),3);

filter_size = 5;
hfs = (filter_size-1)/2;     % half filter size
tos  = icdf('normal',0.5 + 0.9/2, 0, 1)  % ~1.645(fs=11, fs=5) Times of sigma corresponding to cummulative prob = X, so that 2(X-0.5)=0.9.  This is 1-d, does not seem good
tos  = icdf('normal',0.5 + sqrt(0.9/4), 0, 1)  % ~1.949(fs=11, fs=5)  Times of sigma corresponding to cummulative prob = X, so that 4*(x-0.5)^2=0.9  This is 2-d square, also does not work

s_dp = hfs/tos     % when hFilterSize correspond to tos, sigma value in unit of # of data points
h = fspecial('gaussian',filter_size,s_dp)    % the filter with the calculated sigma, and with 90% cummulative prob within area of filterSize.  Filter is square

h = [
    0.0124    0.0263    0.0339    0.0263    0.0124
    0.0263    0.0560    0.0720    0.0560    0.0263
    0.0339    0.0720    0.0925    0.0720    0.0339
    0.0263    0.0560    0.0720    0.0560    0.0263
    0.0124    0.0263    0.0339    0.0263    0.0124]         % this is copied from fitted coefficients
exxA_small= filter2(h,exx0_small);

myplot(exx);caxis([-0.12 0.02])     % (1) raw out put from Vic2D
myplot(exx0_small);caxis([-0.12 0.02])    % (2) 4-point average when reconstructing strain from Vic2D's displacement
myplot(exxA_small);caxis([-0.12 0.02])    % (3) 4-point averaged, then filtered strain, reconstructed from Vic2D's displacement.
myplot(abs(exx-exx0_small) < 1e-3)        %
myplot(abs(exx-exxA_small) < 1e-3)        % Reconstruction from (3) is very similar to direct output (1) from Vic2D












