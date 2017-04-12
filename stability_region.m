a0 = linspace(-4, 4, 20);
b1 = linspace(-4, 4, 20);

w = 1;
h = 1.5;
x = linspace(-4, 4, 100);
y1 = - x * w / (sin(w * h));
y2 = - x / h;
y3 = - w * w / (1 - cos(w * h));
plot(x, y1, 'k', x, y2, 'k', x, y3, '.k');

% for 
%     if(algorithm_g_const(a0, b1, N) > 0) 
%         plot(a0, b0, '+r');
%     else 
%         plot(a0, b0, '.b'); 
%     end
% end