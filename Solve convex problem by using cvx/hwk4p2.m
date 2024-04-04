function [rho,M]= hwk4p2 (Z, a, rhomax, mgiven)
N=520;
cvx_begin
    variable rho(N) 
    variable c(2) 
    variable gammma
    maximize (gammma)
    subject to 
    [a*Z*diag(rho)*Z'-gammma*eye(2) c; c' inv_pos(mgiven)]==semidefinite(3);
    a*sum(rho) == mgiven;
    0 <= rho;
    rho <= rhomax;
    c == (a/mgiven)*Z*rho;
cvx_end
M = a*Z*diag(rho)*Z'-mgiven*c*c'
end
