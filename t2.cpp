//べき乗法による最大固有値の計算
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;

const int N = 6;           // 行列サイズ(書き替え)
const int M = 100;         // 最大反復回数(書き替え)
const double TOL = 1e-8;   // 許容誤差（相対誤差） (書き替え)

//変更不要
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
//変更不要
int main() {
    int mode, digits;

    cout << "出力モードを選んでください:\n";
    cout << "1: 絶対値最大固有値のみ出力\n";
    cout << "2: 固有ベクトルのみ出力\n";
    cout << "3: 両方出力\n";
    cout << "選択（1〜3）: ";
    cin >> mode;

    if (mode < 1 || mode > 3) {
        cout << "1〜3の数字を入力してください。" << endl;
        return 1;
    }

    cout << "出力の小数点以下の桁数（自然数）を入力してください: ";
    cin >> digits;

    if (digits < 0) {
        cout << "自然数を入力してください。" << endl;
        return 1;
    }

    vector<vector<double>> A(N, vector<double>(N));
    vector<double> x(N, 1.0), y(N, 0.0);
    double lambda_old = 0.0, lambda_new = 0.0;

// 行列 A の定義 (書き替え)
    // A[i][j] = 1 + min(i,j)
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            A[i][j] = 1.0 + min(i, j); //<-- 行列の定義(ここを書き換える)

    //以下変更不要
    int iter = 0;
    for (; iter < M; ++iter) {
        for (int i = 0; i < N; ++i) {
            y[i] = 0.0;
            for (int j = 0; j < N; ++j)
                y[i] += A[i][j] * x[j];
        }

        double norm = *max_element(y.begin(), y.end(), [](double a, double b) {
            return fabs(a) < fabs(b);
        });
        for (int i = 0; i < N; ++i)
            x[i] = y[i] / norm;

        lambda_new = compute_rayleigh(A, x);

        if (fabs(lambda_new - lambda_old) / fabs(lambda_new) < TOL)
            break;

        lambda_old = lambda_new;
    }

    cout << "反復回数: " << iter + 1 << endl;
    cout << scientific << setprecision(digits);

    if (mode == 1 || mode == 3)
        cout << "最大固有値 ≈ " << lambda_new << endl;

    if (mode == 2 || mode == 3) {
        cout << "対応する固有ベクトル:" << endl;
        for (int i = 0; i < N; ++i)
            cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}
