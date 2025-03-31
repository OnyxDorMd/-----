import math
from functools import reduce

def factorize(n):
    """Разложение числа на простые множители"""
    print(f"\nФакторизация числа {n}:")
    factors = {}
    # Обработка делителей 2
    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1 #Увеличиваем степень двойки
        n = n // 2
    # Обработка нечетных делителей
    i = 3
    max_factor = math.sqrt(n) + 1 #Верхняя граница для делителей
    while i <= max_factor:
        while n % i == 0:
            factors[i] = factors.get(i, 0) + 1 #Увеличиваем степень i
            n = n // i
            max_factor = math.sqrt(n) + 1
        i += 2 # Переходим к следующему нечетному числу
    if n > 1: # Если остался простой делитель
        factors[n] = factors.get(n, 0) + 1
    print(f"{factors}")
    return factors

def extended_gcd(a, b):
    """Расширенный алгоритм Евклида"""
    if a == 0: # Базовый случай рекурсии
        return (b, 0, 1)
    else:
        g, y, x = extended_gcd(b % a, a)
        res = (g, x - (b // a) * y, y) # Обновляем коэффициенты
        return res

def modinv(a, m):
    """Обратный элемент по модулю"""
    g, x, y = extended_gcd(a, m) #Используем расширенный алгоритм Евклида
    if g != 1: # Обратный элемент существует только, если НОД = 1
        raise Exception('Обратный элемент не существует')
    else:
        res = x % m
        return res

#объединяет два сравнения в одно
def crt_pair(c1, c2):
        m1, r1 = c1
        m2, r2 = c2
        print(f"\nОбъединение сравнений:")
        print(f"x ≡ {r1} mod {m1} и x ≡ {r2} mod {m2}")
        g, p, q = extended_gcd(m1, m2) #Находим коэффициенты х, у
        if (r2 - r1) % g != 0: # Проверка на совместимость
            print("Система не имеет решения!")
            return None
        lcm = m1 // g * m2 # Ищем НОК
        x = (r1 + (r2 - r1) // g * p % (m2 // g) * m1) % lcm # Само решение
        print(f"Общее решение: x ≡ {x} mod {lcm}")
        return (lcm, x) # НОК и решение

def crt(remainders, moduli):
    """Китайская теорема об остатках"""
    print(f"\nКитайская теорема об остатках для:")
    print(f"Остатки: {remainders}")
    print(f"Модули: {moduli}")
        
    result = reduce(crt_pair, zip(moduli, remainders))[1] #Последовательно объединяем все решения
    print(f"\nФинальный результат CRT: x ≡ {result} mod {reduce(lambda x,y: x*y, moduli)}")
    return result

def silver_pohlig_hellman(a, b, q):
    """Алгоритм Сильвера-Полига-Хеллмана"""
    print(f"\n\nВЫЧИСЛЕНИЕ ДИСКРЕТНОГО ЛОГАРИФМА")
    print(f"Находим {a}^x ≡ {b} mod {q}")
    
    # Шаг 1: Факторизация q-1
    factors = factorize(q - 1)
    
    congruences = [] # Список для хранения сравнений
    moduli = [] # Список для хранения модулей
    
    # Шаг 2: Решение подзадач
    for p, alpha in factors.items():
        print(f"\nОбрабатываем простой делитель p={p}, alpha={alpha}")
        
        # Построение таблицы
        print(f"\nПостроение таблицы для p={p}:")
        table = {}
        exponent = (q - 1) // p
        
        for j in range(p):
            val = pow(a, j * exponent, q)
            table[val] = j
            print(f"j={j}: {a}^({j}*{exponent}) mod {q} = {val}")
        
        # Нахождение x mod p^alpha
        x_p = 0
        current_b = b
        
        for k in range(alpha):
            exponent_k = (q - 1) // (p ** (k + 1))
            term = pow(current_b, exponent_k, q) #Вычисляем текущее значение b
            print(f"({current_b})^{exponent_k} mod {q} = {term}")
            
            # Поиск в таблице
            j = table.get(term)
            if j is None:
                raise ValueError(f"Не удалось найти логарифм для p={p}, alpha={alpha}")
            print(f"Найдено в таблице: r{p},{j} = {term} → x{k} = {j}")
            
            # Обновление x_p (частичное решение для текущего р)
            x_p += j * (p ** k)
            
            # Обновление текущего_b
            if k < alpha - 1: # смотрим последний шаг для данного р или нет
                inv = modinv(a, q) # Обратный элемент для а
                exponent_update = j * (p ** k) #текущая часть х
                current_b = (current_b * pow(inv, exponent_update, q)) % q
                
        congruences.append(x_p) #Добавляем в список сравнений
        moduli.append(p ** alpha) #Добавляем в список модулей
        print(f"\nФинальный результат для p={p}: x ≡ {x_p} mod {p**alpha}")
    
    # Шаг 3: КТоО
    print("\nШАГ 3: Китайская теорема об остатках")
    x = crt(congruences, moduli) # Объединение результатов
    
    print(f"\nРЕЗУЛЬТАТ")
    print(f"log_{a}({b}) mod {q} = {x}")
    print(f"Проверка: {a}^{x} mod {q} = {pow(a, x, q)} (ожидалось {b})")
    return x

a = 2
b = 62
q = 181
x = silver_pohlig_hellman(a, b, q)
