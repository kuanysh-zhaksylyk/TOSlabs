

l = 0.1;
a = -20; 
b = 5;

exp_amount = 10;
array_of_Xs = a:l:b; 
n = length(array_of_Xs);
y = zeros(1, n);
array_of_Ys_max = zeros(1, n);
max_x = 0; %#ok<NASGU> 
max_y = -99999; 
index_of_max_y = 0; 
for i_exp = 1:exp_amount
    for i = 1:n
        y(i) = tosmodel5(array_of_Xs(i));
    end
    
    [curr_max_y, index_of_curr_max_y] = max(y);
    if max_y < curr_max_y
        max_y = curr_max_y;
        index_of_max_y = index_of_curr_max_y;
        array_of_Ys_max = y;
    end
end
max_x = array_of_Xs(index_of_max_y); 
fprintf('Max in x: %3.3f \n', max_x);
fprintf('Max in y: %4.4f \n', max_y);
min_y_axis = 15;
max_y_axis = 70; 
plot(array_of_Xs, array_of_Ys_max);
%axis([a b min_y_axis max_y_axis]);
grid on;


% sigma=0.5 
x_init = a:l:b;
n = length(x_init);
y_init = zeros(1, n);
for i = 1:n
    y_init(i) = tosmodel5_new_sigma(x_init(i));
end
[max_y, index_of_max_y] = max(y_init);
max_x = x_init(index_of_max_y);
fprintf('Sigma=0.5\n');
fprintf('Max in x: %1.4f\n', max_x);
fprintf('Max in y: %3.4f\n', max_y);
figure
plot(x_init, y_init);
axis([a b min_y_axis max_y_axis]);
hold on
grid on


steps = 75;
x_nesymm = zeros(1, steps);
y_nesymm = zeros(1, steps);
x_k = -12; 
x_nesymm(1)=x_k
for k = 1:steps
    alpha = 1/(k^(1/3)); 
    gamma = 1/(k); 
    y_sum = tosmodel5_new_sigma(x_k + alpha);
    y_k = tosmodel5_new_sigma(x_k);
    x_k = x_k + (gamma/alpha) * (y_sum - y_k);
    x_nesymm(k+1) = x_k;
    y_nesymm(k+1) = tosmodel5_new_sigma(x_k);
end
x_nesymm
fprintf('Несимметричный алгоритм Кифера–Вольфовица\n');
fprintf('Max in x: %3.4f\n', x_nesymm(steps));
fprintf('Max in y: %4.4f\n', y_nesymm(steps));
plot(x_nesymm, y_nesymm, 'r.');
title('Несимметричный алгоритм Кифера-Вольфица');
axis([a b min_y_axis max_y_axis]);
grid on

% шаги
steps = 150;
x_symm = zeros(1, steps);
y_symm = zeros(1, steps);
x_k = -12;
x_symm(1)=x_k
y_symm(1)=tosmodel5_new_sigma(x_k)
for k = 1:steps
    alpha = 1/(k^(1/3)); 
    gamma = 1/(k); 
    y_sum = tosmodel5_new_sigma(x_k + alpha);
    y_k = tosmodel5_new_sigma(x_k - alpha);
    x_k = x_k + (gamma/alpha) * (y_sum - y_k);
    x_symm(k+1) = x_k;
    y_symm(k+1) = tosmodel5_new_sigma(x_k);
end
x_symm
fprintf('Cимметричный алгоритм Кифера–Вольфовица\n');
fprintf('Max in x: %3.4f\n', x_symm(steps));
fprintf('Max in y: %4.4f\n', y_symm(steps));
hold off
figure
hold on
plot(x_init, y_init);
plot(x_symm, y_symm, 'r.');
title('Симметричный алгоритм Кифера-Вольфица');
%axis([a b min_y_axis max_y_axis]);
grid on

