% This is part of a strain calculation code, not finished yet.
% chenzhe, 2017-09-27

[X,Y] = meshgrid(0:10,0:10);
sigma = rand(size(X));
x = X + rand(11);
y = Y + rand(11);
subset_size = 2;    % half-subset-size
fit_val = 2;
[nR,nC] = size(X);
%%
iR = 3;
iC = 3;

Xcell = X(max([1 iR-subset_size]):min([nR iR+subset_size]), max([1 iC-subset_size]):min([nC iC+subset_size]));
Ycell = Y(max([1 iR-subset_size]):min([nR iR+subset_size]), max([1 iC-subset_size]):min([nC iC+subset_size]));
xcell = x(max([1 iR-subset_size]):min([nR iR+subset_size]), max([1 iC-subset_size]):min([nC iC+subset_size]));
ycell = y(max([1 iR-subset_size]):min([nR iR+subset_size]), max([1 iC-subset_size]):min([nC iC+subset_size]));
scell = sigma(max([1 iR-subset_size]):min([nR iR+subset_size]), max([1 iC-subset_size]):min([nC iC+subset_size]));

X_point = X(iR,iC);
Y_point = Y(iR,iC);

% Determine fit coeffiecients for x,y in terms of X,Y


OK = (scell(:) ~= -1) & (~isnan(xcell(:))) & (~isnan(ycell(:)));

switch fit_val
    case 1
        min_pts = 6;
    case 2
        min_pts = 12;
    case 3
        min_pts = 20;
end

if sum(OK) >= min_pts
    
    XX = Xcell(OK);
    YY = Ycell(OK);
    xx = xcell(OK);
    yy = ycell(OK);
    
    switch fit_val
        case 1
            %  poly11 model
            % (x,y) = aX + bY + c
            A = [XX YY ones(size(YY))];
        case 2
            %  poly22 model
            % (x,y) = aX^2 + bY^2 + cXY + dX + eY + f
            A = [XX.^2 YY.^2 XX.*YY XX YY ones(size(YY))];
        case 3
            %  poly33 model
            % (x,y) = aX^3 + bY^3 + cX^2Y + dXY^2 + eX^2 + fY^2 +
            %              gXY + hX + iY + j
            A = [XX.^3 YY.^3 XX.^2.*YY XX.*YY.^2 XX.^2 ...
                YY.^2 XX.*YY XX YY ones(size(YY))];
    end
    coeff_x = A\xx;
    coeff_y = A\yy;
    switch fit_val
        case 1
            dx_dX = coeff_x(1);
            dx_dY = coeff_x(2);
            dy_dX = coeff_y(1);
            dy_dY = coeff_y(2);
        case 2
            dx_dX = ...
                2*coeff_x(1)*X_point + coeff_x(3)*Y_point + coeff_x(4);
            dx_dY = ...
                2*coeff_x(2)*Y_point + coeff_x(3)*X_point + coeff_x(5);
            dy_dX = ...
                2*coeff_y(1)*X_point + coeff_y(3)*Y_point + coeff_y(4);
            dy_dY = ...
                2*coeff_y(2)*Y_point + coeff_y(3)*X_point + coeff_y(5);
        case 3
            dx_dX = ...
                3*coeff_x(1)*X_point^2 + 2*coeff_x(3)*X_point*Y_point ...
                + coeff_x(4)*Y_point^2 +2*coeff_x(5)*X_point ...
                + coeff_x(7)*Y_point + coeff_x(8);
            dx_dY = ...
                3*coeff_x(2)*Y_point^2 + coeff_x(3)*X_point^2 ...
                + 2*coeff_x(4)*X_point*Y_point +2*coeff_x(6)*Y_point ...
                + coeff_x(7)*X_point + coeff_x(9);
            dy_dX = ...
                3*coeff_y(1)*X_point^2 + 2*coeff_y(3)*X_point*Y_point ...
                + coeff_y(4)*Y_point^2 +2*coeff_y(5)*X_point ...
                + coeff_y(7)*Y_point + coeff_y(8);
            dy_dY = ...
                3*coeff_y(2)*Y_point^2 + coeff_y(3)*X_point^2 ...
                + 2*coeff_y(4)*X_point*Y_point +2*coeff_y(6)*Y_point ...
                + coeff_y(7)*X_point + coeff_y(9);
    end
    
    F = [dx_dX dx_dY; dy_dX dy_dY];
    
end