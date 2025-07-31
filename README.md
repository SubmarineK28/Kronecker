В файле main.cpp содержатся две функции main.   
Первая (закомментированная) с помощью функций createH0(N), createM0(N), createFullHamiltonian(N, beta) полностью строит гамильтониан. Далее функции (см. папку exact_H): 
// double Gauss_tr = Gaussian_trace_estimator(H, k, M);
// double Ray_tr = Rayleigh_trace_estimator(H, k, M);
// double Hutch_tr = Hutchinson_trace_estimator(H, k, M);
// double Un_vec = Unit_vector_estimator(H, k, M);
оценивают след матрицы с помощью различных стохастических методов из статьи,
а функция
// double exact_tr = exact_trace(H, k);
использует точное вычисление следа (от матрицы в степени) с помощью библиотеки Armadillo.
Такой подход работает максимум для 9–10 частиц.

Вторая функция main предназначена для вычисления следа без хранения всей матрицы. Она использует функцию apply_H0, реализующую действие оператора
