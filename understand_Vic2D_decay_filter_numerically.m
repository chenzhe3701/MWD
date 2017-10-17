%
% chenzhe, 2017-10-16 revised
% (1) This script finds out the numeric values of the filters that Vic2D used.
% (2) Two files are required. 
% (3) File-1 should contain the 'exx' processed normally.
% (4) Then, a decay filter of 'filter_size' is applied to 'exx' and it is
% exported as file-2.
% (5) Then, these two files are used to generate new file, and fit filter.

filter_size = 15;
a = load(['WE43_T6_C1_s5_r0c0_fs_',num2str(filter_size),'_a'],'exx');
b = load(['WE43_T6_C1_s5_r0c0_fs_',num2str(filter_size),'_b'],'exx','sigma');
exx1=a.exx;
exx2=b.exx;
sigma=b.sigma;
save(['WE43_T6_C1_s5_r0c0_exx_fs_',num2str(filter_size)],'exx1','exx2','sigma','filter_size')

%%
load(['WE43_T6_C1_s5_r0c0_exx_fs_',num2str(filter_size)],'exx1','exx2','sigma','filter_size');
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