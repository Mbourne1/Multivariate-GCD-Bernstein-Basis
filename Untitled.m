syms x y;

f = (x-2)^2 * (x+y+1)^3 * (y-1);

f_dash = diff(f,x)
expand(f_dash)

f2 = (x-2)* (x+y+1)^2 * (y-1);
f2_dash = diff(f2,x)