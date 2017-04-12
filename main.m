a_min = -4;
a_max = 4;

b_min = -4;
b_max = 4;
n_points = 20;

h = 1.5;
w = 1;

NS = [5, 10, 15, 20, 25];

for N = NS
    x = linspace(a_min, a_max, 1000);
    y1 = - x * w / (sin(w * h));
    y2 = - x / h;
    y3 = - w * w / (1 - cos(w * h));
    plot(x, y1, 'k', x, y2, 'k', x, y3, '.k');
    hold on;
    
    save_dirpath = ['data/', num2str(N), '/'];
    
    mkdir(save_dirpath);
    
    save_fileindex = 0;
    
    for a0 = linspace(a_min, a_max, n_points)
        for b1 = linspace(b_min, b_max, n_points)
            fprintf('-------------------------\n');
            fprintf('### N = %d, a0 = %f, b1 = %f, file=%d\n', N, a0, b1, save_fileindex);
            fprintf('-------------------------\n');
            
            save_filepath = [save_dirpath, num2str(save_fileindex), '/'];
            
            V_min = algorithm_g_const(N, h, w, a0, b1, save_filepath);
            
            if(V_min > 0) 
                plot(a0, b0, '+r');
            else 
                plot(a0, b0, '.b'); 
            end
            hold on;
        end
    end
    
    print([save_dirpath, 'plot'], '-depsc');
    hold off;
end