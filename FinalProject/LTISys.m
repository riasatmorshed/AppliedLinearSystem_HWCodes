function [zeta, omega] = LTISys(a,b)
syms x y
expr1 = -x * y;
expr2 = y*sqrt(1-x^2);

result = solve(expr1 == a, expr2 == b, x,y);
zeta = double(result.x);
omega = double(result.y);
end