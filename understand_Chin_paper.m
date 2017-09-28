% some scripts to help understand the Chin paper
% chenzhe, 2017-09-27

m = [0 1 0]';
n = [1 0 0]';
g = 0.2;

F = eye(3) + g*m*n';
dx_dX = F;
du_dX = F - eye(3);
e_small = 1/2*(du_dX + du_dX');
Vf_Vi = det(F);



P0 = [1,0,0]';   % col vector
P = P0/norm(P0);  % unit vector giving the initial direction of the line
% method-1 basic
lambdaP_square = 0;
for ii=1:3
   for jj=1:3
      for kk=1:3
         lambdaP_square = lambdaP_square +  F(ii,jj)*F(ii,kk)*P(jj)*P(kk);
      end
   end
end
% method-2 matrix
lambdaP_square = (P'*F')*(F*P);     % Eq[5]
lambdaP = lambdaP_square^(1/2);     % ratio of final length to initial length
p = 1/lambdaP * (F*P);              % Eq[6] unit vector giving the final direction of the line



Finv = inv(F);
dX_dx = Finv;

Q0 = [1 0 0]';
Q = Q0/norm(Q0);    % unit vector of the plane normal
% method-1 basic
fQ_square = 0;
for ii=1:3
   for jj=1:3
      for kk=1:3
         fQ_square = fQ_square +  Finv(ii,jj)*Finv(kk,jj)*Q(ii)*Q(kk);
      end
   end
end
% method-2 matrix
fQ_square = Q' * (Finv * Finv') * Q;    % Eq[7]
fQ = fQ_square^(1/2);   % ratio of initial to final perpendicular distance between materials plane, plane unit normal is Q
q = 1/fQ * Finv' * Q;   % Eq[8] unit vector of the final plane normal


%% example

% diagnoal, theoretical expression, but matlab does not calculate this way
% possibly it works better for multiple slip systems
m1 = [0 1 0]';
n1 = [1 0 0]';
m2 = [1 0 0]';
n2 = [0 1 0]';
beta = 2;   % beta = b/a;
gamma = 0.2;
alpha = gamma;

F1 = m1*n1' + beta * m2*n2';
F2 = beta*m2*n2'*m1*n1';
[V,D] = eig(F1);
V * diag(exp(diag(alpha*D))) /(V)     % (1) this sometimes fails
expm(alpha*F1)                    % (2) this is a function that works

F = eye(3) + alpha*F1 + alpha^2*F2 % (3) let (one step shear) a = alpha (total shear)

FA = eye(3) + alpha*m1*n1'
FB = eye(3) + beta*alpha*m2*n2'
F = FB*FA                           % (4) is the same as (3), both assume total shear occurs in one step

% after eq [30] 
F1 = [-2/sqrt(6) 0 0; 0 0 0; 0 -1/sqrt(3) 2/sqrt(6)];
[V,D] = eig(F1)


