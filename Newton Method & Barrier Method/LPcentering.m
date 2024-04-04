function [xc, valc] = LPcentering(A,b,c,x_0)
ALPHA = 0.01;
BETA = 0.5;
EPSILON = 1e-6;
MAXITERS = 100;
if (min(x_0) <= 0) || (norm(A*x_0 - b) > 1e-3) 
    fprintf('FAILED');
    nu_star = []; xc = []; lambda_hist=[];
    return;
end
m = length(b);
n = length(x_0);
x = x_0; lambda_hist = [];valc=[];
for iter = 1:MAXITERS
    val = c'* x-sum(log(x))
    valc = [valc, val];
    H = diag(x.^(-2));
    g = c - x.^(-1);
    w = (A*diag(x.^2)*A')\(-A*diag(x.^2)*g);
    dx = -diag(x.^2)*(A'*w + g);
    lambdasqr = -g'*dx; % dx'*H*dx;
    lambda_hist = [lambda_hist lambdasqr/2];
    if lambdasqr/2 <= EPSILON break; end
    t = 1; while min(x+t*dx) <= 0 t = BETA*t; end
    while c'*(t*dx)-sum(log(x+t*dx))+sum(log(x))-ALPHA*t*g'*dx> 0
          t = BETA*t;end
      x = x + t*dx;
end
if iter == MAXITERS
      fprintf('ERROR: MAXITERS reached.\n');
      xc = []; nu_star = [];
  else
      xc = x;
      nu_star = w;
end
semilogy(lambda_hist,'bo-')
xlabel('iters')
ylabel('lambdasqr/2')