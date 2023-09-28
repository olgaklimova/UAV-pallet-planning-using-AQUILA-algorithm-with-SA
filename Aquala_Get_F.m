%function [LB,UB,Dim,F_obj] = Aquala_Get_F(xc, yc, zc, r, Li_x, Li_y, Li_z, Linext_x, Linext_y, Linext_z)
function [LB,UB,Dim,F_obj] = Aquala_Get_F(xc, yc, r, Li_x, Li_y, Linext_x, Linext_y)
F_obj = @F1;
LB=-350;
UB=350;
Dim = 2;
    
% F1
function o = F1(x)

%центр препятствия
%c_pr = [xc yc zc];
c_pr = [xc yc];
R = r;
Ds = 0.9; %значение данного параметра в статье 0.906
Us = 0.3; %значение данного параметра в статье 0.302


% %направляющий вектор прямой (отрезка LiLi+1)
% %s = [Linext_x - Li_x; Linext_y - Li_y; Linext_z - Li_z];
% %s = [Linext_x - Li_x; Linext_y - Li_y];
% % s = [x(1) - Li_x; x(2) - Li_y; x(3) - Li_z];
% s = [x(1) - Li_x; x(2) - Li_y];
% %точка лежащая на прямой
% M1 = [Li_x Li_y Li_z];
% %M1 = [Li_x Li_y];
% %расстояние от c_pr до M1
% c_prM1 = c_pr - M1;
% %векторное произведение c_prM1 на s
% c_prM1s = cross(c_prM1*s);
% %нахождение расстояния
% %D = abs(c_prM1s)/abs(sqrt(s(1).^2 + s(2).^2 + s(3).^2));
% D = abs(c_prM1s)/abs(sqrt(s(1).^2 + s(2).^2));

%расстояние от точки до прямой
A = (x(2)-Li_y)/(x(1)-Li_x);
B = -1;
C = Li_y-Li_x*(x(2)-Li_y)/(x(1)-Li_x);
D = abs(A*c_pr(1)+B*c_pr(2)+C)/sqrt(A^2+B^2);

%

%значение 
% параментр для целевой функции
if D > R + Ds
o_pr = 0;
elseif D < R + Us
o_pr = inf;
else  
o_pr = Ds - D;
end

% целевая функция состоящая из минимизации расстояния и возможности огибания препятствий
%o= sqrt((Li_x - x(1)).^2 + (Li_y - x(2)).^2 + (Li_z - x(3)).^2) + sqrt((Linext_x - x(1)).^2 + (Linext_y - x(2)).^2 + (Linext_z - x(3)).^2) + o_pr;
%o= sqrt((Li_x - x(1)).^2 + (Li_y - x(2)).^2 + (Li_z - x(3)).^2) + sqrt((Linext_x - x(1)).^2 + (Linext_y - x(2)).^2 + (Linext_z - x(3)).^2);
if sqrt((c_pr(1)-x(1))^2+(c_pr(2)-x(2))^2)> R + Ds
    o= sqrt((Li_x - x(1)).^2 + (Li_y - x(2)).^2) + sqrt((Linext_x - x(1)).^2 + (Linext_y - x(2)).^2);
else
    o= sqrt((Li_x - x(1)).^2 + (Li_y - x(2)).^2) + sqrt((Linext_x - x(1)).^2 + (Linext_y - x(2)).^2) + o_pr;
%o= sqrt((Li_x - x(1)).^2 + (Li_y - x(2)).^2) + sqrt((Linext_x - x(1)).^2 + (Linext_y - x(2)).^2);
end
end

% целевая функция учитываючая только расстояние
% function o = F1(x)
% o= sqrt((2 - x(1)).^2 + (2 - x(2)).^2) + sqrt((5 - x(1)).^2 + (3 - x(2)).^2);
end