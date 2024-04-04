function [Y] = hw4q3 (W)
n=40;
cvx_begin sdp
    variable X(n,n) symmetric
    minimize(-trace(W*(ones(n,n)-X))/4)
    subject to
    diag(X)==ones(n,1)
    X>=0
    X == semidefinite(n);
cvx_end
Y=X
end


