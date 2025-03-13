#include <iostream>
#include <ceres/ceres.h>
#include <glog/logging.h>

// 定義一個簡單的代價函數
struct CostFunctor {
    template <typename T>
    bool operator()(const T* const x, T* residual) const {
        // f(x) = 10 - x，我們要最小化 (f(x))^2
        residual[0] = 10.0 - x[0];
        return true;
    }
};

int main(int argc, char** argv) {
    // 初始化 glog
    google::InitGoogleLogging(argv[0]);

    // 初始參數值
    double x = 0.5;
    double initial_x = x;

    // 建立問題
    ceres::Problem problem;

    // 使用 auto-diff 創建代價函數
    ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);

    // 添加殘差項到問題中
    problem.AddResidualBlock(cost_function, nullptr, &x);

    // 設定求解器選項
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;

    // 解決問題
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // 輸出結果
    std::cout << summary.BriefReport() << "\n";
    std::cout << "初始值 x = " << initial_x << "\n";
    std::cout << "最終值 x = " << x << "\n";
    
    // 理想結果應該是 x 約等於 10
    std::cout << "預期值 x = 10" << std::endl;

    return 0;
}