//任意の個数(1000以下)の与えられたミューオンの寿命データ(ヘッダーがない単純な形式のCSVファイル)を元に、寿命の指数分布を仮定し、尤度関数を定義し、その後、最尤推定量を計算してミューオンの寿命を推定するROOTライブラリ使ったC++のプログラム

#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

// データをファイルから読み込む関数
std::vector<double> read_data(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Unable to open file");
    }

    std::vector<double> data;
    double value;
    while (file >> value) {
        data.push_back(value);
    }

    if (data.size() > 1000) {
        throw std::runtime_error("Data size exceeds the limit of 1000");
    }

    return data;
}

// 対数尤度関数を計算する関数
double neg_log_likelihood(const double* params) {
    double mu = params[0];
    double sum = 0.0;
    for (double value : data) {
        sum += log(mu) - mu * value;
    }
    return -sum;
}

int main() {
    try {
        std::vector<double> data = read_data("muon_lifetimes.csv");

        // 最小化オブジェクトを作成
        ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

        // 最大関数呼び出し回数、最大反復回数、収束許容範囲を設定
        minimizer->SetMaxFunctionCalls(10000);
        minimizer->SetMaxIterations(1000);
        minimizer->SetTolerance(0.001);

        // 対数尤度関数を最小化するためのFunctorを設定
        ROOT::Math::Functor f(&neg_log_likelihood, 1);
        minimizer->SetFunction(f);

        double step[1] = {0.01};
        double variable[1] = {1.0};  // 初期値

        // 最小化する変数を設定
        minimizer->SetLimitedVariable(0, "mu", variable[0], step[0], 0, 0);

        // 最小化を実行
        minimizer->Minimize();

        // 推定された平均寿命を取得して表示
        const double *xs = minimizer->X();
        std::cout << "Estimated mean lifetime: " << xs[0] << std::endl;

        delete minimizer;
    } catch (std::runtime_error& error) {
        std::cerr << "Error: " << error.what() << std::endl;
    }

    return 0;
}
