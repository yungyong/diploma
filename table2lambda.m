function f = table2lambda(X, Y)
    pp = spline(X, Y);
    f = @(x) ppval(pp, x);
end

