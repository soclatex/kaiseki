// jacobi法を用いて、n次の線形方程式系を解くプログラム
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;

const int N = 6;                // 行列サイズ(書き替え)
const int M = 100;              // 最大反復回数(書き替え)
const double TOL = 1e-8;        // 許容誤差（相対誤差） (書き替え)

int main(){
    int digits;
    cout << "出力の小数点以下の桁数（自然数）を入力してください: ";
    cin >> digits;

    if (digits < 0) {
        cout << "自然数を入力してください。" << endl;
        return 1;
    }

    vector<vector<double>> A(N, vector<double>(N, 0.0));
    vector<double> b(N, 1.0);
    vector<double> x_old(N, 1.0); // 初期値 x^(0)，すべての成分を1.0に設定
    vector<double> x_new(N, 0.0);

    // 行列 A の定義 (書き替え)
    for (int i = 0; i < N; ++i) {
        A[i][i] = 5.0; // 対角成分を5.0に設定
        if (i > 0) A[i][i - 1] = 1.0;
        if (i < N - 1) A[i][i + 1] = 1.0; // 対角の隣接要素を1.0に設定

        if (i > 1) A[i][i - 2] = 1.0;
        if (i < N - 2) A[i][i + 2] = 1.0; // 対角の2つ離れた(隣接)要素を1.0に設定
    }
//以下変更不要
    int iter = 0;
    for (; iter < M ; ++iter) {
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j)
                if (j != i)
                    sum += A[i][j] * x_old[j];
            x_new[i] = (b[i] - sum) / A[i][i];
        }

        double max_diff = 0.0, max_val = 0.0;
        for (int i = 0; i < N; ++i) {
            max_diff = max(max_diff, fabs(x_new[i] - x_old[i]));
            max_val = max(max_val, fabs(x_new[i]));
        }

        if (max_diff / max_val < TOL)
            break;

        x_old = x_new;
    }

    cout << "反復回数: " << iter + 1 << endl;
    cout << scientific << setprecision(digits);
    cout << "近似解の第1成分 x1 = " << x_new[0] << endl;

    return 0;
}
