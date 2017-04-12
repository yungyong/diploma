function [xa, ya] = matr_Lyap(a0, a1, b1, h, w)

% data prepare 
L = [a0, a1, b1; -a1, -a0, -b1; 1, -1, 0];
e = expm(h*L);
    syms z(t) v(t) x1(t) ;
    z_h = e(1,1)*z(0) + e(1,2)*v(0) + e(1,3)*x1(0);
    v_h = e(2,1)*z(0) + e(2,2)*v(0) + e(2,3)*x1(0);
    x1_h = e(3,1)*z(0) + e(3,2)*v(0) + e(3,3)*x1(0);
    
    c1 = a0 + a0*e(2,1) + a1*e(1,1) + b1*e(3,1);
    c2 = a1 + a0*e(2,2) + a1*e(1,2) + b1*e(3,2);
    c3 = b1 + a0*e(2,3) + a1*e(1,3) + b1*e(3,3);
    
    A = [1 - e(2,1), -e(2,2), - e(2,3); -e(3,1), -e(3,2), 1 - e(3,3); c1, c2, c3];
    b = [0; 0; -w];
    C = rref([A b]);
    y1 = C(1:3,4:4);

%     %ode 
%        eqs = [
%         diff(z, t) == z*a0 + v*a1 + x1*b1 ,...
%         diff(v, t) == -a1*z - a0*v - b1*x1,...
%         diff(x1, t) == z - v
%         ];
% 
%     % initial conditions
%        cond = [
%         z(0) == y1(1),...
%         x1(0) == y1(3),...
%         v(0) == y1(2)
%         ];
% 
%     s = dsolve(eqs, cond);
%     z = simplify(s.z);
    
% system solve 
fun = @(x,y) [y(1)*a0 + y(2)*a1 + y(3)*b1; -a1*y(1) - a0*y(2) - b1*y(3); y(1) - y(2)];
[xa, ya] = ode45(fun, [0 h], [y1(1) y1(2) y1(3)]);
ya = ya(:, 1);

%%%
xr = xa(2:end);
xr = flip(-xr);

yr = ya(2:end);
yr = flip(yr);

xa = [xr; xa]';
ya = [yr; ya]';
end
