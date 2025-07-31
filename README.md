В файле main.cpp содержатся две функции main.   <br>
Первая (закомментированная) с помощью функций createH0(N), createM0(N), createFullHamiltonian(N, beta) полностью строит гамильтониан. Далее функции (см. папку exact_H): <br>
// double Gauss_tr = Gaussian_trace_estimator(H, k, M); <br>
// double Ray_tr = Rayleigh_trace_estimator(H, k, M); <br>
// double Hutch_tr = Hutchinson_trace_estimator(H, k, M); <br>
// double Un_vec = Unit_vector_estimator(H, k, M); <br>
оценивают след матрицы с помощью различных стохастических методов из статьи,  <br>
а функция <br>
// double exact_tr = exact_trace(H, k); <br>
использует точное вычисление следа (от матрицы в степени) с помощью библиотеки Armadillo. <br>
Такой подход работает максимум для 9–10 частиц. <br>
 <br>
Вторая функция main предназначена для вычисления следа без хранения всей матрицы. Она использует функцию apply_H0, реализующую действие оператора <br>
H0 на вектор z (без явного построения матрицы). Она была предложена чатом GPT, и использует методы квантового программирования, которые я еще не совсем понимаю.
