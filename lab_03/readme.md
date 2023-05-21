## Результат

Пока так 

<img src="img/img_14.png" width=500px>

## Краткая попытка собрать все, что мы знаем

### Физический смысл задачи

Есть труба, по которой течет горячий воздух. Вокруг трубы воздух холодный. Нам нужно определить температуру стенки трубки. 

<img src="img/img_01.svg" width=500px>

Как выглядит стенка в разрезе (цветом показано уменьшение температуры от теплого к холодному):

<img src="img/img_02.svg" width=500px>

### Формульный смысл

```
Лирическое отступление на тему квазилинейных уравнений.

Квазилинейное уравнение - это уравнение, в котором нелинейные 
члены зависят от произведений производных неизвестной функции.

Нелинейные члены - члены в какой-то не первой степени.
```

Нам дано квазилинейное уравнение от `T(r)`.

<img src="img/img_03.svg" width=600px>

Мы хотим построить функцию `T(r)` [хотя в итоге на графике будет `T(z)`]. То есть хотим понять, как изменяется температура в зависимости от удаления.

Для этого у нас есть два краевых условия.

<img src="img/img_04.svg" width=500px>

📍 **Зачем нужны краевые условия?**

``` Они 
Краевые условия позволяют узнать начальные прогоночные коэффициенты и yn
для обратной прогонки.
```

📍 **Что делать, если краевые условия не линейные?**

```
Сейчас краевые условия линейны (так как все в первой степени).
Если бы они были нелинейными, то есть три варианта решения:
	1. Метод простых итераций
	2. Метод Ньютона
	3. Изменение направления прогонки
Подробнее в приложении
```

### Еще раз итог дано = что нужно сделать

Нужно посторить график зависимости T от r.

### Способ решения задачи

Важно понимать, что мы решаем задачу следующего вида и отталкиваемся от всех этих `u''`, `p` и `f`.

<img src="img/img_05.svg" width=500px>

Используем **разностную схему** - заменим производные на их разностные аналоги.

<img src="img/img_06.svg" width=350px>

Обращаем внимание на `i` - то есть мы можем заменить только конкретные значения -> получаем набор уравнений. Если у нас есть n точек, то неизвестных игреков n + 1 (так как во второй производной используется i + 1). Вопрос, где нам найти n + 1 уравнение, чтобы это дело получить.

n - 1 уравнение у нас уже есть в тот момент, когда мы использовали разностную схему. Добавляем два уравнения - краевых условия и получаем n + 1. Ура.

<img src="img/img_07.svg" width=500px>

Эта система является трехдиагональной матрицей, которая имеет следующий вид:

<img src="img/img_08.svg" width=200px>

Эту матрицу можно описать следующей системой уравнений:

<img src="img/img_09.svg" width=500px>

Будем считать, что все коэффициенты здесь мы знаем (`A`, `B` ...). И мы их действительно можем выразить из наших начальных условий. Вопрос который нужно решить на этом этапе -- как зная все коэффициенты, найти y. 

### Метод прогонки

```
Метод прогонки - это численный метод решения трехдиагональных систем линейных алгебраических уравнений (СЛАУ), которые могут быть записаны в виде матричного уравнения Ax = b, где A - трехдиагональная матрица коэффициентов, x - вектор неизвестных и b - вектор правой части.
```

Будем считать, что между y есть следующая зависимость.

<img src="img/img_12.svg" width=250px>

Как найти хси и ету 

<img src="img/img_10.svg" width=500px>

Как найти первые кси и ету

<img src="img/img_11.svg" width=500px>

Отсюда мы можем выразить yn и используя основную прогоночную формулу (см выше) посчитать y назад. 

Теперь мы знаем один вариант y, но как понять, что он правильный?

### Метод простых итераций

```
Метод простой итерации - это численный метод решения систем линейных алгебраических уравнений (СЛАУ), которые могут быть записаны в виде матричного уравнения Ax = b, где A - матрица коэффициентов, x - вектор неизвестных и b - вектор правой части.
```

После рассчета уточним ими коэффициенты А, B и так далее и посчитаем y заново. Будем продолжать, пока разница не станет меньше некоторой ошибки. Все...

<img src="img/img_13.svg" width=500px>


