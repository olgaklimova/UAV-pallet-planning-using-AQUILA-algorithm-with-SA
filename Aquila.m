clc;
clear;
close all;


% Построение графика зависимости значения целевой функции от количества препятствий
% x = [566.52 686.73 687.72];
% y = [5 8 10];
% plot(y, x)

%MAP = Aquala_Map([8 8], [inf inf inf inf 4]);
figure(1)
hold on;

%ЗНАЧЕНИЯ ПУТИ И ПРЕПЯТСТВИЙ ИЗ СТАТЬИ
%генерация коньена размером из статьи
% map1 = Aquala_Map([3500 3500], [0 3500 3500 0 3500]);
map1 = Aquala_Map([350 350], [0 350 350 0 350]);

max_h = max(map1(:));
min_h = min(map1(:));
hh = max_h+abs(min_h);
z = [hh hh];

% ТЕСТ № 1
% %координаты пути
% % x = [1000 3000];
% % y = [1000 2000];
% x = [100 300];
% y = [100 200];
% surf(map1,'edgecolor','none')
% plot3(x,y,z,'k-*');
% 
% %препятствия (5 штук)
% k = 5;
% %координаты центров
% % xc = [1400 2500 2200 1700 2000];
% % yc = [1300 2000 1200 2000 1600];
% xc = [140 250 220 170 200];
% yc = [130 200 120 200 160];
% zc = [hh hh hh hh hh hh];
% %радиусы
% r = [10 15 10 20 7.5];

% % ТЕСТ № 2
% %координаты пути
% x = [50 250];
% y = [80 250];
% surf(map1,'edgecolor','none')
% plot3(x,y,z,'k-*');
% 
% %препятствия (8 штук)
% k = 8;
% %координаты центров
% xc = [150 250 220 170 200 100 80 220];
% yc = [130 200 120 200 160 150 100 230];
% zc = [hh hh hh hh hh hh hh hh];
% %радиусы
% r = [10 15 10 20 7.5 10 8 7.5];

% ТЕСТ № 3
%координаты пути
x = [50 250];
y = [80 250];
surf(map1,'edgecolor','none')
plot3(x,y,z,'k-*');

%препятствия (10 штук)
k = 10;
%координаты центров
xc = [150 250 220 170 200 100 80 220 140 215];
yc = [130 200 120 200 160 150 100 230 160 190];
zc = [hh hh hh hh hh hh hh hh hh hh];
%радиусы
r = [10 15 10 20 15 10 8 7.5 8 7.5];

%построение препятствий
for i = 1:k
plot3(xc(i),yc(i),zc(i),'b-*');
g = map1(xc(i),yc(i));
t = g:pi/100:max_h+abs(min_h)*pi;
xt1 = r(i)*sin(5*t) + xc(i);
yt1 = r(i)*cos(5*t) + yc(i);
plot3(xt1,yt1,t)
end

% %РАНДОМНЫЕ ЗНАЧЕНИЯ ПУТИ И ПРЕПЯТСТВИЙ (для одного препятствия)
% %обычная генерация карты
% % map1 = Aquala_Map([100 100], [inf inf inf inf 100]);
% 
% max_h = max(map1(:));
% min_h = min(map1(:));
% x = randi(100,2,1);
% y = randi(100,2,1);
% z = zeros(size(x));
% for i = 1:length(x)
%     %z(i) = map1(y(i),x(i)) + (max_h + abs(min_h));
%     %самое простое задание высоты
%     z(i) = max_h;
%     %иное определение высоты точек
%     %if (map1(y(i),x(i)) < 0)
%     %z(i) = map1(y(i),x(i)) + 1.5*max_h
%     %else
%     %z(i) = map1(y(i),x(i)) + max_h;
%     %end
%     %траектория под определенным углом (высота подъема точек тоже рандомная - но иногда пересекает ландшафт)
%     %z(i) = map1(y(i),x(i))+(max_h-map1(y(i),x(i)))*rand([1,1]);
% end
% surf(map1,'edgecolor','none')
% plot3(x,y,z,'k-*');
% 
% %формируем препятствие
% 
% %определяем наиболшее, наименьшее значение в векторе x
% maxx = max(x(1),x(2));
% minx = min(x(1),x(2));
% 
% %берем рандомное значение координаты x на прямой для центра препятствий
% xc = randi([minx,maxx],1,1);
% % нахождение y из уравнения прямой
% yc = ((xc - x(1))*(y(2)-y(1)))/(x(2)-x(1)) + y(1);
% % нахождение z из уравнения прямой
% zc = ((yc - y(1))*(z(2)-z(1)))/(y(2)-y(1)) + z(1);
% 
% while sqrt((xc - x(1)).^2 + (yc - y(1)).^2 + (zc - z(1)).^2) <= 1 || sqrt((xc - x(2)).^2 + (yc - y(2)).^2 + (zc - z(2)).^2) <= 1
%     %берем рандомное значение координаты x на прямой для центра препятствий
%     xc = randi([minx,maxx],1,1);
%     % нахождение y из уравнения прямой
%     yc = ((xc - x(1))*(y(2)-y(1)))/(x(2)-x(1)) + y(1);
%     % нахождение z из уравнения прямой
%     zc = ((yc - y(1))*(z(2)-z(1)))/(y(2)-y(1)) + z(1);
% end
% 
% plot3(xc,yc,zc,'b-*');
% 
% % определяем радиус от 1 до 5
% r = randi([1,5],1,1);
% %r = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2+(z(1)-z(2))^2)/10/2;
% while r > sqrt((xc - x(1)).^2 + (yc - y(1)).^2 + (zc - z(1)).^2) || r > sqrt((xc - x(2)).^2 + (yc - y(2)).^2 + (zc - z(2)).^2)
%     r = randi([1,5],1,1);
% end

% g = map1(round(yc),round(xc));
% t = g:pi/100:max_h+abs(min_h)*pi;
% xt1 = r*sin(5*t) + xc;
% yt1 = r*cos(5*t) + yc;
% plot3(xt1,yt1,t)

x = [50 100];
y = [80 90];
% Разбиение пути на n+1 отрезок n точками
n = 30; %в статье данный параметр равен 30
Li_x = zeros(1,n);
Li_y = zeros(1,n);
Li_z = zeros(1,n);

for i = 1:n
lamb = i/(n+1-i);
Li_x(i) = (x(1) + lamb*x(2))/(1+lamb);
Li_y(i) = (y(1) + lamb*y(2))/(1+lamb);
Li_z(i) = (z(1) + lamb*z(2))/(1+lamb);

%plot3(Li_x(i),Li_y(i),Li_z(i),'r-o');
end

%Применение Aquala optimization algorithm

Solution_no=20;  % N особей в алгоритме Aquala
M_Iter=1000;     % Количество итераций: T в алгоритме Aquala

Best_P_x = zeros(1,n+1);
Best_P_y = zeros(1,n+1);
Best_P_z = zeros(1,n+1);

%[LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, zc, r, x(1), y(1), z(1), Li_x(1), Li_y(1), Li_z(1)); 
[LB,UB,Dim,F_obj]= Aquala_Get_F(xc(1), yc(1), r(1), x(1), y(1), Li_x(1), Li_y(1)); 
[~,Best_P,~] = Aquila_AO(Solution_no,M_Iter,LB,UB,Dim,F_obj); 
Best_P_x(1) = Best_P(1);
Best_P_y(1) = Best_P(2);
Best_P_z(1) = z(1);
%plot3(Best_P_x(1),Best_P_y(1),Best_P_z(1),'m-p')

for i = 1:n-1
%[LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, zc, r, Li_x(i), Li_y(i), Li_z(i), Li_x(i+1), Li_y(i+1), Li_z(i+1)); 
[LB,UB,Dim,F_obj]= Aquala_Get_F(xc(1), yc(1), r(1), Li_x(i), Li_y(i), Li_x(i+1), Li_y(i+1)); 
[Best_FF(i),Best_P,conv] = Aquila_AO(Solution_no,M_Iter,LB,UB,Dim,F_obj); 
Best_P_x(i+1) = Best_P(1);
Best_P_y(i+1) = Best_P(2);
Best_P_z(i+1) = Li_z(i);
% if(Best_P_x(i+1) == 0)
%  sqrt((Li_x - x(1)).^2 + (Li_y - x(2)).^2) + sqrt((Linext_x - x(1)).^2 + (Linext_y - x(2)).^2))
% end
%plot3(Best_P_x(i+1),Best_P_y(i+1),Li_z(i+1),'m-p')
end

%[LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, zc, r, Li_x(n), Li_y(n), Li_z(n), x(2), y(2), z(2)); 
[LB,UB,Dim,F_obj]= Aquala_Get_F(xc(1), yc(1), r(1), Li_x(n), Li_y(n), x(2), y(2)); 
[Best_FF(n+1),Best_P,conv] = Aquila_AO(Solution_no,M_Iter,LB,UB,Dim,F_obj); 
Best_P_x(n+1) = Best_P(1);
Best_P_y(n+1) = Best_P(2);
Best_P_z(n+1) = z(2);
%plot3(Best_P_x(n+1),Best_P_y(n+1),Best_P_z(n+1),'m-p')

plot3(Best_P_x,Best_P_y,Best_P_z,'m-p')

display('The optimal path is built');
% display(['The best-obtained solution by AO is : ', num2str(Best_P)]);

Best_Path = 0;
for i = 1:n
Best_Path = Best_Path + (sqrt((Best_P_x(i+1) - Best_P_x(i)).^2 + (Best_P_y(i+1) - Best_P_y(i)).^2 + (Best_P_z(i+1) - Best_P_z(i)).^2));
end
display(['The length of path found by AO is : ', num2str(Best_Path)]);
display(['The best optimal value of the objective funciton found by AO is : ', num2str(sum(Best_FF))]);

x = [100 250];
y = [90 250];
% Разбиение пути на n+1 отрезок n точками
n = 100; %в статье данный параметр равен 30
Li_x = zeros(1,n);
Li_y = zeros(1,n);
Li_z = zeros(1,n);

for i = 1:n
lamb = i/(n+1-i);
Li_x(i) = (x(1) + lamb*x(2))/(1+lamb);
Li_y(i) = (y(1) + lamb*y(2))/(1+lamb);
Li_z(i) = (z(1) + lamb*z(2))/(1+lamb);

%plot3(Li_x(i),Li_y(i),Li_z(i),'r-o');
end

%Применение Aquala optimization algorithm

Solution_no=20;  % N особей в алгоритме Aquala
M_Iter=1000;     % Количество итераций: T в алгоритме Aquala

Best_P_x = zeros(1,n+1);
Best_P_y = zeros(1,n+1);
Best_P_z = zeros(1,n+1);

%[LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, zc, r, x(1), y(1), z(1), Li_x(1), Li_y(1), Li_z(1)); 
[LB,UB,Dim,F_obj]= Aquala_Get_F(xc(1), yc(1), r(1), x(1), y(1), Li_x(1), Li_y(1)); 
[~,Best_P,~] = Aquila_AO(Solution_no,M_Iter,LB,UB,Dim,F_obj); 
Best_P_x(1) = Best_P(1);
Best_P_y(1) = Best_P(2);
Best_P_z(1) = z(1);
%plot3(Best_P_x(1),Best_P_y(1),Best_P_z(1),'m-p')

for i = 1:n-1
%[LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, zc, r, Li_x(i), Li_y(i), Li_z(i), Li_x(i+1), Li_y(i+1), Li_z(i+1)); 
[LB,UB,Dim,F_obj]= Aquala_Get_F(xc(1), yc(1), r(1), Li_x(i), Li_y(i), Li_x(i+1), Li_y(i+1)); 
[Best_FF(i),Best_P,conv] = Aquila_AO(Solution_no,M_Iter,LB,UB,Dim,F_obj); 
Best_P_x(i+1) = Best_P(1);
Best_P_y(i+1) = Best_P(2);
Best_P_z(i+1) = Li_z(i);
% if(Best_P_x(i+1) == 0)
%  sqrt((Li_x - x(1)).^2 + (Li_y - x(2)).^2) + sqrt((Linext_x - x(1)).^2 + (Linext_y - x(2)).^2))
% end
%plot3(Best_P_x(i+1),Best_P_y(i+1),Li_z(i+1),'m-p')
end

%[LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, zc, r, Li_x(n), Li_y(n), Li_z(n), x(2), y(2), z(2)); 
[LB,UB,Dim,F_obj]= Aquala_Get_F(xc(1), yc(1), r(1), Li_x(n), Li_y(n), x(2), y(2)); 
[Best_FF(n+1),Best_P,conv] = Aquila_AO(Solution_no,M_Iter,LB,UB,Dim,F_obj); 
Best_P_x(n+1) = Best_P(1);
Best_P_y(n+1) = Best_P(2);
Best_P_z(n+1) = z(2);
%plot3(Best_P_x(n+1),Best_P_y(n+1),Best_P_z(n+1),'m-p')

plot3(Best_P_x,Best_P_y,Best_P_z,'m-p')

display('The optimal path is built');
% display(['The best-obtained solution by AO is : ', num2str(Best_P)]);

Best_Path = 0;
for i = 1:n
Best_Path = Best_Path + (sqrt((Best_P_x(i+1) - Best_P_x(i)).^2 + (Best_P_y(i+1) - Best_P_y(i)).^2 + (Best_P_z(i+1) - Best_P_z(i)).^2));
end
display(['The length of path found by AO is : ', num2str(Best_Path)]);
display(['The best optimal value of the objective funciton found by AO is : ', num2str(sum(Best_FF))]);

% %поправка точек, которые лежат внутри препятствия
% %%
% %syms x_var y_var;
% for i = 1:n
% if sqrt((xc(1)-Li_x(i))^2+(yc-Li_y(i))^2)<= r
%     b = r;
%     c = sqrt((xc-Li_x(i))^2+(yc-Li_y(i))^2);
%     a = sqrt(b^2-c^2);
% %     A = atan(a/b);
% %     B = atan((y(1)-y(2))/(x(1)-x(2)));
% %     Li_x(i) = xc+cos(rad2deg(A+B))*b;
% %     Li_y(i) = yc+sin(rad2deg(A+B))*b;
%     x1 =(1/2)*((yc-Li_y(i))*sqrt(-(-xc^2+2*Li_x(i)*xc-Li_x(i)^2+(-r+a-yc+Li_y(i))*(-r+a+yc-Li_y(i)))*(-xc^2+2*Li_x(i)*xc-Li_x(i)^2+(r+a-yc+Li_y(i))*(r+a+yc-Li_y(i)))*(xc-Li_x(i))^2)+(xc^3-xc^2*Li_x(i)+(Li_y(i)^2-2*yc*Li_y(i)-r^2+yc^2+a^2-Li_x(i)^2)*xc-Li_x(i)*(a^2-r^2-Li_x(i)^2-Li_y(i)^2+2*yc*Li_y(i)-yc^2))*(xc-Li_x(i)))/((xc-Li_x(i))*(xc^2-2*Li_x(i)*xc+Li_x(i)^2+(yc-Li_y(i))^2));
%     
%     y1 = (-sqrt(-(-xc^2+2*Li_x(i)*xc-Li_x(i)^2+(-r+a-yc+Li_y(i))*(-r+a+yc-Li_y(i)))*(-xc^2+2*Li_x(i)*xc-Li_x(i)^2+(r+a-yc+Li_y(i))*(r+a+yc-Li_y(i)))*(xc-Li_x(i))^2)+yc^3-yc^2*Li_y(i)+(a^2+xc^2-r^2+Li_x(i)^2-2*Li_x(i)*xc-Li_y(i)^2)*yc+Li_y(i)^3+(Li_x(i)^2-2*Li_x(i)*xc+r^2-a^2+xc^2)*Li_y(i))/(2*yc^2-4*yc*Li_y(i)+2*Li_y(i)^2+2*(xc-Li_x(i))^2);
%     
%     
%     x2 = (1/2)*((-yc+Li_y(i))*sqrt(-(-xc^2+2*Li_x(i)*xc-Li_x(i)^2+(-r+a-yc+Li_y(i))*(-r+a+yc-Li_y(i)))*(xc-Li_x(i))^2*(-xc^2+2*Li_x(i)*xc-Li_x(i)^2+(r+a-yc+Li_y(i))*(r+a+yc-Li_y(i))))+(xc-Li_x(i))*(xc^3-xc^2*Li_x(i)+(yc^2-2*yc*Li_y(i)+Li_y(i)^2+a^2-r^2-Li_x(i)^2)*xc-Li_x(i)*(-r^2-Li_x(i)^2+a^2-yc^2+2*yc*Li_y(i)-Li_y(i)^2)))/((xc^2-2*Li_x(i)*xc+Li_x(i)^2+(yc-Li_y(i))^2)*(xc-Li_x(i)));
%     
%     y2 = (sqrt(-(xc-Li_x(i))^2*(-xc^2+2*Li_x(i)*xc-Li_x(i)^2+(r+a+yc-Li_y(i))*(r+a-yc+Li_y(i)))*(-xc^2+2*Li_x(i)*xc-Li_x(i)^2+(-r+a+yc-Li_y(i))*(-r+a-yc+Li_y(i))))+yc^3-yc^2*Li_y(i)+(a^2+xc^2-r^2+Li_x(i)^2-2*Li_x(i)*xc-Li_y(i)^2)*yc+Li_y(i)^3+(Li_x(i)^2-2*Li_x(i)*xc+r^2-a^2+xc^2)*Li_y(i))/(2*yc^2-4*yc*Li_y(i)+2*Li_y(i)^2+2*(xc-Li_x(i))^2);
%     Li_x(i) = x1;
%     Li_y(i) = y1;
% end
% %plot3(Li_x(i),Li_y(i),Li_z(i),'b-o');
% end
% end

%%
% % Применение Aquala optimization algorithm
% 
% Solution_no=20;  % N особей в алгоритме Aquala
% M_Iter=1000;     % Количество итераций: T в алгоритме Aquala
% 
% Best_P_x = zeros(1,n+1);
% Best_P_y = zeros(1,n+1);
% Best_P_z = zeros(1,n+1);
% 
% %[LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, zc, r, x(1), y(1), z(1), Li_x(1), Li_y(1), Li_z(1)); 
% [LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, r, x(1), y(1), Li_x(1), Li_y(1)); 
% [~,Best_P,~] = Aquila_AO(Solution_no,M_Iter,LB,UB,Dim,F_obj); 
% Best_P_x(1) = Best_P(1);
% Best_P_y(1) = Best_P(2);
% Best_P_z(1) = z(1);
% %plot3(Best_P_x(1),Best_P_y(1),Best_P_z(1),'m-p')
% 
% for i = 1:n-1
% %[LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, zc, r, Li_x(i), Li_y(i), Li_z(i), Li_x(i+1), Li_y(i+1), Li_z(i+1)); 
% [LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, r, Li_x(i), Li_y(i), Li_x(i+1), Li_y(i+1)); 
% [Best_FF(i),Best_P,conv] = Aquila_AO(Solution_no,M_Iter,LB,UB,Dim,F_obj); 
% Best_P_x(i+1) = Best_P(1);
% Best_P_y(i+1) = Best_P(2);
% Best_P_z(i+1) = Li_z(i);
% % if(Best_P_x(i+1) == 0)
% %  sqrt((Li_x - x(1)).^2 + (Li_y - x(2)).^2) + sqrt((Linext_x - x(1)).^2 + (Linext_y - x(2)).^2))
% % end
% %plot3(Best_P_x(i+1),Best_P_y(i+1),Li_z(i+1),'m-p')
% end
% 
% %[LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, zc, r, Li_x(n), Li_y(n), Li_z(n), x(2), y(2), z(2)); 
% [LB,UB,Dim,F_obj]= Aquala_Get_F(xc, yc, r, Li_x(n), Li_y(n), x(2), y(2)); 
% [Best_FF(n+1),Best_P,conv] = Aquila_AO(Solution_no,M_Iter,LB,UB,Dim,F_obj); 
% Best_P_x(n+1) = Best_P(1);
% Best_P_y(n+1) = Best_P(2);
% Best_P_z(n+1) = z(2);
% %plot3(Best_P_x(n+1),Best_P_y(n+1),Best_P_z(n+1),'m-p')
% 
% plot3(Best_P_x,Best_P_y,Best_P_z,'m-p')
% 
% display('The optimal path is built');
% % display(['The best-obtained solution by AO is : ', num2str(Best_P)]);
% 
% Best_Path = 0;
% for i = 1:n
% Best_Path = Best_Path + (sqrt((Best_P_x(i+1) - Best_P_x(i)).^2 + (Best_P_y(i+1) - Best_P_y(i)).^2 + (Best_P_z(i+1) - Best_P_z(i)).^2));
% end
% display(['The length of path found by AO is : ', num2str(Best_Path)]);
% display(['The best optimal value of the objective funciton found by AO is : ', num2str(sum(Best_FF))]);