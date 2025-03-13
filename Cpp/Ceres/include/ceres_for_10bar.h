// ceres_for_10bar.h
#ifndef CERES_FOR_10BAR_H
#define CERES_FOR_10BAR_H

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include <vector>
#include <cmath>

constexpr double PI = 3.14159265358979323846;
constexpr double E = 200e9;  // Young's modulus
constexpr double LENGTH_UNIT = 9.14;  // Element unit length
constexpr double MAX_DISPLACEMENT = 0.02;  // Maximum allowed displacement
constexpr double MAX_STRESS = 250e6;  // Maximum allowed stress
constexpr double DENSITY = 7860;  // Material density

// Structure to hold truss element data
struct TrussElement {
    int node1;
    int node2;
    double angle;  // in radians
    int areaType;  // 0 for first group, 1 for second group
    double length;
};

// Class to handle FEM computation
class TrussFEMSolver {
public:
    TrussFEMSolver();
    
    // Compute displacement and stress for given radii
    bool Solve(const double r1, const double r2, 
               std::vector<double>& displacements,
               std::vector<double>& stresses);
    
    // Helper functions
    Eigen::MatrixXd ComputeStiffnessMatrix(double angle);
    double ComputeWeight(double r1, double r2) const;

private:
    void InitializeTruss();
    
    // Truss data
    std::vector<TrussElement> elements_;
    std::vector<int> bcNodes_;  // Boundary condition nodes
    std::vector<double> forceVector_;
    int numNodes_;
};

// Cost functor for weight objective
struct WeightCostFunctor {
    WeightCostFunctor(TrussFEMSolver* solver) : solver_(solver) {}

    template <typename T>
    bool operator()(const T* const r, T* residual) const {
        // Compute weight as the objective using templated types for autodiff
        T r1_squared = r[0] * r[0];
        T r2_squared = r[1] * r[1];
        T sqrt2 = T(std::sqrt(2.0));
        
        // Calculate weight: (6*r1^2 + 4*r2^2*sqrt(2))*pi*density*length
        residual[0] = (T(6.0) * r1_squared + T(4.0) * r2_squared * sqrt2) * 
                      T(PI) * T(DENSITY) * T(LENGTH_UNIT);
        return true;
    }

private:
    TrussFEMSolver* solver_;
};

// 使用数值微分的位移约束函数
struct NumericDisplacementConstraintFunctor {
    NumericDisplacementConstraintFunctor(TrussFEMSolver* solver) : solver_(solver) {}

    bool operator()(const double* const r, double* residual) const {
        std::vector<double> displacements;
        std::vector<double> stresses;
        
        if (!solver_->Solve(r[0], r[1], displacements, stresses)) {
            // 求解失败，设置一个大的残差
            residual[0] = 1e10;
            return false;
        }
        
        // 计算最大位移
        double maxDisp = 0.0;
        for (size_t i = 0; i < displacements.size()/2; ++i) {
            double dx = displacements[2*i];
            double dy = displacements[2*i+1];
            double disp = std::sqrt(dx*dx + dy*dy);
            maxDisp = std::max(maxDisp, disp);
        }
        
        // 实际不等式约束: maxDisp <= MAX_DISPLACEMENT
        // 转换为残差: maxDisp - MAX_DISPLACEMENT <= 0
        residual[0] = std::max(0.0, maxDisp - MAX_DISPLACEMENT);
        
        return true;
    }

private:
    TrussFEMSolver* solver_;
};

// 数值微分的应力约束函数
struct NumericStressConstraintFunctor {
    NumericStressConstraintFunctor(TrussFEMSolver* solver, int elementIndex) 
        : solver_(solver), elementIndex_(elementIndex) {}

    bool operator()(const double* const r, double* residual) const {
        std::vector<double> displacements;
        std::vector<double> stresses;
        
        if (!solver_->Solve(r[0], r[1], displacements, stresses)) {
            // 求解失败，设置一个大的残差
            residual[0] = 1e10;
            return false;
        }
        
        // 获取当前元素的应力
        double stress = stresses[elementIndex_];
        
        // 实际不等式约束: |stress| <= MAX_STRESS
        // 转换为残差: |stress| - MAX_STRESS <= 0
        double abs_stress = std::abs(stress);
        residual[0] = std::max(0.0, abs_stress - MAX_STRESS);
        
        // 为了让残差的量级更适合优化器，可以进行缩放
        residual[0] *= 1e-8;  // 应力单位通常较大
        
        return true;
    }

private:
    TrussFEMSolver* solver_;
    int elementIndex_;
};

// 主优化函数
bool OptimizeTruss(double& r1, double& r2);

#endif // CERES_FOR_10BAR_H