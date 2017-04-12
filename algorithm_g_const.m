function final = algorithm_g_const(N, h, w, a0, b1)
a1 = 0;
del = h/N;

[X, Y] = lyapunov_matrix(a0, a1, b1, h, w);
u = table2lambda(X, Y);

L = zeros(1, N);
M = zeros(1, N);

parfor i = 1:N 
    I1_L = integral(@(s) u(-s - del * i) .* a1 .* (1 + s / del), -del, 0);
    I1_M = integral(@(s) u(-s - del * i) .* a1 .* (- s / del), -del, 0);

    I2_internal_scalar = @(s) integral(@(xi) u(xi - s + (N - i) * del) .* b1, -h, s - (N - i) * del);
    I2_internal = @(S) arrayfun(I2_internal_scalar, S);

    I2_L = integral(@(s) I2_internal(s) .* (1 + s / del), -del, 0);
    I2_M = integral(@(s) I2_internal(s) .* (- s / del), -del, 0);

    L(i) = I1_L + I2_L;
    M(i) = I1_M + I2_M;
end

fprintf('L, M are ready...\n');
fprintf('=================\n');

P = zeros(N);
Q = zeros(N);
R = zeros(N);

parfor i = 1:N
    for j = 1:N
        I1_internal = @(s1, s2) a1^2 .* u(s1 - s2 - (i-j) .* del); 
        I1_P = integral(@(s1) integral(@(s2)I1_internal(s1, s2) .* (1 + s2/del), -del, 0) .* (1 + s1/del), -del, 0);
        I1_Q = integral(@(s1) integral(@(s2)I1_internal(s1, s2) .* (- s2/del), -del, 0) .* (1 + s1/del), -del, 0);
        I1_R = integral(@(s1) integral(@(s2)I1_internal(s1, s2) .* (s2/del), -del, 0) .* (s1/del), -del, 0); 
        
        I2_internal_scalar = @(s1, s2) 2 .* integral(@(xi2) a1 .* u(h + s1 - s2 - (i - j) * del - xi2) .* b1, -h, s2 - (N - i) * del);
        I2_internal = @(S1, S2) arrayfun(I2_internal_scalar, S1, S2);
        I2_P = integral(@(s1) integral(@(s2)I2_internal(s1, s2) .* (1 + s2/del), -del, 0) .* (1 + s1/del), -del, 0);
        I2_Q = integral(@(s1) integral(@(s2)I2_internal(s1, s2) .* (- s2/del), -del, 0) .* (1 + s1/del), -del, 0);
        I2_R = integral(@(s1) integral(@(s2)I2_internal(s1, s2) .* (s2/del), -del, 0) .* (s1/del), -del, 0);
        
        I3_insternal_scalar_subinternal_scalar = @(s1, s2, xi2) integral(@(xi1) b1 .* u(xi2 + s1 - s2 - (i - j) .* del - xi1) .* b1,...
            -h, s1 - (N - j) * del);
        I3_insternal_scalar_subinternal = @(s1, s2, XI2) arrayfun(@(xi2) I3_insternal_scalar_subinternal_scalar(s1, s2, xi2), XI2);
        I3_internal_scalar = @(s1, s2) integral(@(xi2) I3_insternal_scalar_subinternal(s1, s2, xi2), -h, s2 - (N - i) * del);
        I3_internal = @(S1, S2) arrayfun(I3_internal_scalar, S1, S2);
        I3_P = integral(@(s1) integral(@(s2) I3_internal(s1, s2) .* (1 + s2/del), -del, 0) .* (1 + s1/del), -del, 0);
        I3_Q = integral(@(s1) integral(@(s2) I3_internal(s1, s2) .* (- s2/del), -del, 0) .* (1 + s1/del), -del, 0);
        I3_R = integral(@(s1) integral(@(s2) I3_internal(s1, s2) .* (s2/del), -del, 0) .* (s1/del), -del, 0);
        
        P(i, j) = I1_P + I2_P + I3_P;
        Q(i, j) = I1_Q + I2_Q + I3_Q;
        R(i, j) = I1_R + I2_R + I3_R;
        
        fprintf('P, Q, R: (%d, %d) is ready...\n', i, j);
    end
end

fprintf('=================\n');
        
c = zeros(1, N);
for i = 1:(N - 1)
    c(i) = L(N - i) + M(N - i + 1) + P(N, N - i) + Q(N, N - i + 1);
end
c(N) = M(1) + Q(N, 1);

fprintf('c is ready...\n');
fprintf('=================\n');

C = zeros(N);
for i = 1:(N - 1)
    for j = 1:(N - 1)
        C(i, j) = P(N - i, N - j) + 2 * Q(N - i, N - j + 1) + R(N - i + 1, N - j + 1);
        C(N, j) = R(1, N - j + 1);
    end
    C(i, N) = 2 * Q(N - 1, 1) + R(N - i + 1, 1);
end

C(N, N) = R(1, 1);

fprintf('C is ready...\n');
fprintf('=================\n');

psi_1 = 0;

parfor i = 1:N
    I1_psi = - integral(@(s) abs(u(-s - del * i)) .* abs(a1) .* (s .^ 2 - del * s) .* (abs(a0) + abs(a1) + abs(b1) * h) .^ 2, -del, 0);

    I2_internal_psi = @(s) integral(@(xi) abs(u(xi - s + (N - i) * del)) .* abs(b1), -h, s - (N - i) * del);
    I2_internal = @(S) arrayfun(I2_internal_psi, S);

    I2_psi = - integral(@(s) I2_internal(s) .* (s .^ 2 - del * s) .* (abs(a0) + abs(a1) + abs(b1) * h) .^ 2, -del, 0);
     
    psi_1 = psi_1 + I1_psi + I2_psi;

    fprintf('PSI1 %d is ready...\n', i);
end

fprintf('=================\n');

psi_2 = 0;
const = 1/2 * (abs(a0) + abs(a1) + abs(b1) * h) ^ 2;
parfor i = 1:N
    for j = 1:N
        I1_internal = @(s1, s2) a1^2 .* abs(u(s1 - s2 - (i-j) .* del)); 
        I1_psi = integral(@(s1) integral(@(s2)I1_internal(s1, s2) .* (2 + const * (s2 .^ 2 - del * s2)), -del, 0) .* (s1 .^ 2 - del * s1), -del, 0);
        
        I2_internal_scalar = @(s1, s2) 2 .* integral(@(xi2) abs(a1) .* abs(u(h + s1 - s2 - (i - j) * del - xi2)) .* abs(b1), -h, s2 - (N - i) * del);
        I2_internal = @(S1, S2) arrayfun(I2_internal_scalar, S1, S2);
        I2_psi = integral(@(s1) integral(@(s2)I2_internal(s1, s2) .* (2 + const * (s2 .^ 2 - del * s2)), -del, 0) .* (s1 .^ 2 - del * s1), -del, 0);
        
        I3_insternal_scalar_subinternal_scalar = @(s1, s2, xi2) integral(@(xi1) abs(b1) .* abs(u(xi2 + s1 - s2 - (i - j) .* del - xi1)) .* abs(b1),...
            -h, s1 - (N - j) * del);
        I3_insternal_scalar_subinternal = @(s1, s2, XI2) arrayfun(@(xi2) I3_insternal_scalar_subinternal_scalar(s1, s2, xi2), XI2);
        I3_internal_scalar = @(s1, s2) integral(@(xi2) I3_insternal_scalar_subinternal(s1, s2, xi2), -h, s2 - (N - i) * del);
        I3_internal = @(S1, S2) arrayfun(I3_internal_scalar, S1, S2);
        I3_psi = integral(@(s1) integral(@(s2) I3_internal(s1, s2) .* (2 + const * (s2 .^ 2 - del * s2)), -del, 0) .* (s1 .^ 2 - del * s1), -del, 0);
        
        psi_2 = psi_2 + I1_psi + I2_psi + I3_psi;   
        
        fprintf('PSI2 (%d, %d) is ready...\n', i, j);
    end
end
psi_2 = psi_2 * (- const);

fprintf('=================\n');

c_scal = u(0) + 2 * L(N) + P(N, N) + psi_1 + psi_2;
bound = diag(eye(N));
x0 = zeros(N, 1);
minfun = @(x) c_scal + 2 * c * x + x' * C * x;
answer = fmincon(minfun, x0, [], [], [], [], -1 .* bound, bound);

final = minfun(answer);
end
