
%{ Laguirre Polynomials as basis functions for use in LSM Monte Carlo. Here maximum k (5) basis functions used. }%

function A = BasisFunctions(X,k)

if k == 1
    A = [ones(size(X)) (1-X)];
elseif k == 2
    A = [ones(size(X)) (1-X) 1/2*(2-4*X+X.^2)];
elseif k == 3
    A = [ones(size(X)) (1-X) 1/2*(2-4*X+X.^2) ...
        1/6*(6-18*X+9*X.^2-X.^3)];
elseif k == 4
    A = [ones(size(X)) (1-X) 1/2*(2-4*X+X.^2) ...
        1/6*(6-18*X+9*X.^2-X.^3) ...
        1/24*(24-96*X+72*X.^2-16*X.^3+X.^4)];
else if k == 5
    A = [ones(size(X)) (1-X) 1/2*(2-4*X+X.^2) ...
            1/6*(6-18*X+9*X.^2-X.^3) ...
            1/24*(24-96*X+72*X.^2-16*X.^3+X.^4) ...
            1/120*(120-600*X+600*X.^2-200*X.^3+25*X.^4-X.^5)];
else
    error('Too many basis functions requested');
    end
end

