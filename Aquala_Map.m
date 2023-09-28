%Генерация карты

function MAP = Aquala_Map(m_size, m_prop)
% Формирование основы карты
% Базовая карта должна быть квадратной после чего обрезается до нужного размера заданного заданием
n=ceil(log2(max(max(m_size(1), m_size(2)))-1));
MAP = inf(2^n+1, 2^n+1);
% Коэффициент разброса вероятности
h = m_prop(5);
% Значения в углах карты
if ismember(Inf, m_prop(1:4)) || ismember(-Inf, m_prop(1:4))
    % Если среди начальных значений есть невозможные (+/- бесконечность)
    % то генерируем случайные значения
    MAP(1,1) = rand;
    MAP(size(MAP,1),1) = rand;
    MAP(1,size(MAP,2)) = rand;
    MAP(size(MAP,1),size(MAP,2)) = rand;
else
    % В противном случае берем заданные значения
    MAP(1,1) = m_prop(1);
    MAP(size(MAP,1),1) = m_prop(2);
    MAP(1,size(MAP,2)) = m_prop(3);
    MAP(size(MAP,1),size(MAP,2)) = m_prop(4);
end

% Генерация остальных значений карты
for i = 1 : n
    % ширина (глубина) шага вычисления значений функции (шаг квадратов)
    m = 2^(n-i);
    for j = 1:2:2^i
        for k = 1:2:2^i
            if MAP(j*m+1, k*m+1) == inf
                MAP(j*m+1, k*m+1) = (MAP((j+1)*m+1, (k+1)*m+1) + MAP((j-1)*m+1, (k-1)*m+1) + MAP((j+1)*m+1, (k-1)*m+1) + MAP((j-1)*m+1, (k+1)*m+1))/4 + sign(rand-0.5)*rand*h;
            end
        end
    end
    
    % Шаг ромбов
    for j = 0 : 2^i
        % Проверка наличия смещения на текущей заполняемой строке 
        % (наличие первного вычисленного значения)
        if mod(j,2)==0
            sk=1;
        else
            sk=0;
        end
        for k = sk : 2 : 2^i
            if MAP(j*m+1, k*m+1) == inf
                s = [];
                p = 0;
                % Проверка границ карты
                if (j+1)*m+1 <= 2^n+1
                    s = [s, MAP((j+1)*m+1, (k)*m+1)];
                    p = p + 1;
                end
                if (j-1)*m+1 >= 1
                    s = [s, MAP((j-1)*m+1, (k)*m+1)];
                    p = p + 1;
                end
                if (k+1)*m+1 <= 2^n+1
                    s = [s, MAP((j)*m+1, (k+1)*m+1)];
                    p = p + 1;
                end
                if (k-1)*m+1 >= 1
                    s = [s, MAP((j)*m+1, (k-1)*m+1)];
                    p = p + 1;
                end
                MAP(j*m+1, k*m+1) = sum(s)/p + sign(rand-0.5)*rand*h;
            end
        end
    end
    % Снижение влияния разброса
    h = h/2;
end

% % Обрезка карты по размеру
if size(MAP,1) ~= m_size(1) && size(MAP,2) ~= m_size(2)
    MAP = MAP(1+round((size(MAP,1)-m_size(1))/2):m_size(1)+round((size(MAP,1)-m_size(1))/2), 1+round((size(MAP,2)-m_size(2))/2):m_size(2)+round((size(MAP,2)-m_size(2))/2));
end

end