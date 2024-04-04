function [x_star, history, value] = LPsolverf(A,b,c,x_0)
T_0 = 1;
MU = 20;
EPSILON = 1e-3; 
n = length(x_0);
t = T_0;
x = x_0;
history = [];
while(1)
    [x_star, valc] = LPcentering(A,b,t*c,x);
    x = x_star;
    gap = n/t;
    history = [history [length(valc); gap]];
    if gap < EPSILON break; end
    t = MU*t;
end
value = c'*x_star;
[xx, yy] = stairs(cumsum(history(1,:)),history(2,:));
semilogy(xx,yy,'bo-');
xlabel('iters')
ylabel('gap')
cvx_begin
    variable x(n)
    minimize(c'*x)
    subject to 
    A*x == b
    x >= 0 
cvx_end

fprintf('\n\nOptimal value found by barrier method:\n');
value
fprintf('Optimal value found by CVX:\n');
cvx_optval
fprintf('Duality gap from barrier method:\n');
gap

