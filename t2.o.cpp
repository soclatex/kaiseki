//絶対最大固有値を求めるためのべき乗法の実装
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;

const int N = 6;              // 行列サイズ
const int M = 100;             // 最大反復回数
const double TOL = 1e-8;     // 許容誤差（相対誤差）

// Rayleigh商 λ ≈ (xᵀAx) / (xᵀx)
double compute_rayleigh(const vector<vector<double>>& A, const vector<double>& x) {
    double num = 0.0, den = 0.0;
    vector<double> Ax(N, 0.0);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            Ax[i] += A[i][j] * x[j];
    for (int i = 0; i < N; ++i) {
        num += x[i] * Ax[i];
        den += x[i] * x[i];
    }
    return num / den;
}

int main() {
    vector<vector<double>> A(N, vector<double>(N));
    vector<double> x(N, 1.0);
    vector<double> y(N, 0.0);
    double lambda_old = 0.0, lambda_new = 0.0;

    // 行列 A の定義
    // A[i][j] = 1 + min(i,j)
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            A[i][j] = 1.0 + min(i, j);  //<-- 行列の定義(ここを書き換える)
            

    int iter = 0;
    for (; iter < M ; ++iter) {
        // y = A * x
        for (int i = 0; i < N; ++i) {
            y[i] = 0.0;
            for (int j = 0; j < N; ++j)
                y[i] += A[i][j] * x[j];
        }

        // 正規化（y -> x）
        double norm = *max_element(y.begin(), y.end(), [](double a, double b) {
            return fabs(a) < fabs(b);
        });
        for (int i = 0; i < N; ++i)
            x[i] = y[i] / norm;

        // 固有値の更新（Rayleigh商）
        lambda_new = compute_rayleigh(A, x);

        if (fabs(lambda_new - lambda_old) / fabs(lambda_new) < TOL)
            break;

        lambda_old = lambda_new;
    }

    // 出力
    cout << "反復回数: " << iter + 1 << endl;
    cout << scientific << setprecision(7);　// 小数点以下7桁で出力(書き替え)
    cout << "最大固有値 ≈ " << lambda_new << endl;

    cout << "対応する固有ベクトル:" << endl;
    for (int i = 0; i < N; ++i)
        cout << "x[" << i << "] = " << x[i] << endl;

    return 0;
}
