import java.util.Scanner;
import java.util.Locale;
import java.util.function.Function;

public class MESProjekt {

    //parametry calkowania numerycznego (Gauss-Legendre 2-punktowy)
    private static final double[] ROOTS = {-0.5773502691896257, 0.5773502691896257}; //punkty
    private static final double[] WEIGHTS = {1.0, 1.0};                             //wagi

    public static void main(String[] args) {

        //konfiguracja
        Scanner scanner = new Scanner(System.in);
        System.out.println("podaj liczbe elementow n: ");
        int n;

        //obsluga argumentu uruchomieniowego
        if (args.length > 0) {
            try {
                n = Integer.parseInt(args[0]);
            } catch (NumberFormatException e) {
                n = scanner.nextInt();
            }
        } else {
            n = scanner.nextInt();
        }

        if (n < 2) {
            System.out.println("n musi byc >= 2");
            return;
        }

        double L_len = 2.0;
        double h = L_len / n;
        int nodes = n + 1;

        double[][] K = new double[nodes][nodes];
        double[] F = new double[nodes];

        System.out.println("obliczenia numeryczne dla n=" + n);

        //agregacja z calkowaniem numerycznym
        for (int i = 0; i < n; i++) {
            double x_start = i * h;
            double x_end = (i + 1) * h;

            //dla kazdego elementu licze 4 calki (macierz 2x2)
            // K_local_00 (dN1 * dN1)
            double k00 = integrateGauss(x_start, x_end, x -> k(x) * dN1(h) * dN1(h));
            // K_local_01 (dN1 * dN2)
            double k01 = integrateGauss(x_start, x_end, x -> k(x) * dN1(h) * dN2(h));
            // K_local_10 (dN2 * dN1)
            double k10 = integrateGauss(x_start, x_end, x -> k(x) * dN2(h) * dN1(h));
            // K_local_11 (dN2 * dN2)
            double k11 = integrateGauss(x_start, x_end, x -> k(x) * dN2(h) * dN2(h));

            //agregacja do macierzy globalnej
            K[i][i]     += k00;
            K[i][i+1]   += k01;
            K[i+1][i]   += k10;
            K[i+1][i+1] += k11;
        }

        //warunki brzegowe

        //lewy brzeg (Robin): u'(0) + u(0) = 20
        //modyfikacja wynikajaca ze sformulowania wariacyjnego
        K[0][0] -= 1.0;  //czlon -u(0)v(0)
        F[0]    -= 20.0; //czlon -20v(0)

        //prawy brzeg (Dirichlet): u(2) = 0
        //zeruje ostatni wiersz i ustawiam 1 na przekatnej
        for (int j = 0; j < nodes; j++) {
            K[n][j] = 0.0;
        }
        K[n][n] = 1.0;
        F[n] = 0.0; //wymuszona wartosc 0

        //rozwiazanie
        double[] u = gaussianElimination(K, F);

        //wyniki
        System.out.println("\nskopiuj do excela:");
        for (int i = 0; i < nodes; i++) {
            System.out.printf("%.4f\t%.4f\n", (i * h), u[i]);
        }
    }

    //funkcje pomocnicze (fizyka i ksztalt)

    //parametr materialowy k(x)
    private static double k(double x) {
        if (x <= 1.0) return 1.0;
        else return 2.0;
    }

    //pochodna pierwszej funkcji ksztaltu (stala = -1/h)
    private static double dN1(double h) {
        return -1.0 / h;
    }

    //pochodna drugiej funkcji ksztaltu (stala = 1/h)
    private static double dN2(double h) {
        return 1.0 / h;
    }

    //silnik calkowania numerycznego
    //implementacja kwadratury Gaussa-Legendre (dwa punkty)
    private static double integrateGauss(double a, double b, Function<Double, Double> func) {
        double center = (b + a) / 2.0;
        double halfLength = (b - a) / 2.0;
        double sum = 0.0;

        for (int i = 0; i < ROOTS.length; i++) {
            //transformacja z przedzialu [-1,1] na [a,b]
            double x_real = center + halfLength * ROOTS[i];
            double weight = WEIGHTS[i];

            //sumowanie: waga * wartosc funkcji * jakobian (halfLength)
            sum += weight * func.apply(x_real);
        }
        return sum * halfLength;
    }

    //eliminacja Gaussa
    public static double[] gaussianElimination(double[][] A, double[] b) {
        int N = b.length;
        double[] x = new double[N];
        for (int k = 0; k < N - 1; k++) {
            for (int i = k + 1; i < N; i++) {
                double factor = A[i][k] / A[k][k];
                b[i] -= factor * b[k];
                for (int j = k; j < N; j++) {
                    A[i][j] -= factor * A[k][j];
                }
            }
        }
        for (int i = N - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }
}
