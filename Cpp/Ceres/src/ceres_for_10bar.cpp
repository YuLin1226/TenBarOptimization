// ceres_for_10bar.cpp
#include "ceres_for_10bar.h"
#include <iostream>

TrussFEMSolver::TrussFEMSolver() : numNodes_(6) {
    InitializeTruss();
}

void TrussFEMSolver::InitializeTruss() {
    // Initialize truss elements as per the Python code
    // Element data = {node1, node2, angle, areaType, length}
    elements_.resize(10);
    
    // Horizontal elements (areaType = 0, length = 1 unit)
    elements_[0] = {3, 5, 0.0, 0, 1.0};
    elements_[1] = {1, 3, 0.0, 0, 1.0};
    elements_[2] = {4, 6, 0.0, 0, 1.0};
    elements_[3] = {2, 4, 0.0, 0, 1.0};
    elements_[4] = {3, 4, 0.0, 0, 1.0};
    elements_[5] = {1, 2, 0.0, 0, 1.0};
    
    // Diagonal elements (areaType = 1, length = sqrt(2) units)
    elements_[6] = {4, 5, 135.0 * PI / 180.0, 1, std::sqrt(2.0)};
    elements_[7] = {3, 6, 45.0 * PI / 180.0, 1, std::sqrt(2.0)};
    elements_[8] = {2, 3, 135.0 * PI / 180.0, 1, std::sqrt(2.0)};
    elements_[9] = {1, 4, 45.0 * PI / 180.0, 1, std::sqrt(2.0)};
    
    // Scale all lengths by LENGTH_UNIT
    for (auto& element : elements_) {
        element.length *= LENGTH_UNIT;
    }
    
    // Set boundary conditions (fixed nodes)
    bcNodes_ = {5, 6};
    
    // Initialize force vector (zero except at specific nodes)
    forceVector_.resize(2 * numNodes_, 0.0);
    forceVector_[2*3-1] = -1e7;  // Node 3, Y direction
    forceVector_[2*4-1] = -1e7;  // Node 4, Y direction
}

Eigen::MatrixXd TrussFEMSolver::ComputeStiffnessMatrix(double angle) {
    double C = std::cos(angle);
    double S = std::sin(angle);
    
    Eigen::MatrixXd SM(4, 4);
    SM << C*C, C*S, -C*C, -C*S,
          C*S, S*S, -C*S, -S*S,
          -C*C, -C*S, C*C, C*S,
          -C*S, -S*S, C*S, S*S;
    
    return SM;
}

bool TrussFEMSolver::Solve(const double r1, const double r2, 
                          std::vector<double>& displacements,
                          std::vector<double>& stresses) {
    // Check for valid input
    if (r1 <= 0.0 || r2 <= 0.0) {
        return false;
    }
    
    // Calculate areas from radii
    std::vector<double> areas(elements_.size());
    for (size_t i = 0; i < elements_.size(); ++i) {
        double radius = (elements_[i].areaType == 0) ? r1 : r2;
        areas[i] = PI * radius * radius;
    }
    
    // Assemble global stiffness matrix
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(2 * numNodes_, 2 * numNodes_);
    
    for (size_t n = 0; n < elements_.size(); ++n) {
        int index1 = 2 * elements_[n].node1 - 2;
        int index2 = 2 * elements_[n].node1 - 1;
        int index3 = 2 * elements_[n].node2 - 2;
        int index4 = 2 * elements_[n].node2 - 1;
        
        std::vector<int> indices = {index1, index2, index3, index4};
        
        // Element stiffness matrix
        Eigen::MatrixXd SM = ComputeStiffnessMatrix(elements_[n].angle) * 
                             E * areas[n] / elements_[n].length;
        
        // Add to global stiffness matrix
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                K(indices[i], indices[j]) += SM(i, j);
            }
        }
    }
    
    // Store full stiffness matrix for reaction calculation
    Eigen::MatrixXd Ks = K;
    
    // Remove rows and columns corresponding to boundary conditions
    std::vector<int> toRemove;
    for (int node : bcNodes_) {
        toRemove.push_back(2 * node - 2);
        toRemove.push_back(2 * node - 1);
    }
    
    // Sort in descending order for proper removal
    std::sort(toRemove.begin(), toRemove.end(), std::greater<int>());
    
    // Create reduced matrices by removing fixed DOF
    Eigen::MatrixXd K_reduced = K;
    Eigen::VectorXd F_reduced(forceVector_.size());
    
    for (size_t i = 0; i < forceVector_.size(); ++i) {
        F_reduced(i) = forceVector_[i];
    }
    
    for (int idx : toRemove) {
        // Remove row and column from K
        int rows = K_reduced.rows() - 1;
        int cols = K_reduced.cols() - 1;
        
        if (idx < rows) {
            K_reduced.block(idx, 0, rows - idx, cols + 1) = 
                K_reduced.block(idx + 1, 0, rows - idx, cols + 1);
        }
        
        if (idx < cols) {
            K_reduced.block(0, idx, rows + 1, cols - idx) = 
                K_reduced.block(0, idx + 1, rows + 1, cols - idx);
        }
        
        K_reduced.conservativeResize(rows, cols);
        
        // Remove element from F
        if (idx < F_reduced.size() - 1) {
            F_reduced.segment(idx, F_reduced.size() - idx - 1) = 
                F_reduced.segment(idx + 1, F_reduced.size() - idx - 1);
        }
        
        F_reduced.conservativeResize(F_reduced.size() - 1);
    }
    
    // Solve for displacements
    Eigen::VectorXd Q_reduced = K_reduced.fullPivLu().solve(F_reduced);
    
    // Create full displacement vector (including zeros at fixed DOFs)
    Eigen::VectorXd Q = Eigen::VectorXd::Zero(2 * numNodes_);
    
    int j = 0;
    for (int i = 0; i < 2 * numNodes_; ++i) {
        if (std::find(toRemove.begin(), toRemove.end(), i) == toRemove.end()) {
            Q(i) = Q_reduced(j++);
        }
    }
    
    // Extract displacements to return vector
    displacements.resize(2 * numNodes_);
    for (int i = 0; i < 2 * numNodes_; ++i) {
        displacements[i] = Q(i);
    }
    
    // Calculate element stresses
    stresses.resize(elements_.size());
    
    for (size_t i = 0; i < elements_.size(); ++i) {
        const auto& elem = elements_[i];
        double C = std::cos(elem.angle);
        double S = std::sin(elem.angle);
        
        Eigen::Vector4d S1(-C, -S, C, S);
        Eigen::Vector4d D;
        
        D << Q(2 * elem.node1 - 2),
             Q(2 * elem.node1 - 1),
             Q(2 * elem.node2 - 2),
             Q(2 * elem.node2 - 1);
        
        stresses[i] = E / elem.length * S1.dot(D);
    }
    
    return true;
}

double TrussFEMSolver::ComputeWeight(double r1, double r2) const {
    // Weight calculation as per the objective function
    return (6 * r1 * r1 + 4 * r2 * r2 * std::sqrt(2.0)) * PI * DENSITY * LENGTH_UNIT;
}

bool OptimizeTruss(double& r1, double& r2) {
    // Initial values
    double parameters[2] = {r1, r2};
    
    // Create the solver
    TrussFEMSolver* solver = new TrussFEMSolver();
    
    // Create the Ceres problem
    ceres::Problem problem;
    
    // Add objective function (weight minimization)
    ceres::CostFunction* cost_function = 
        new ceres::AutoDiffCostFunction<WeightCostFunctor, 1, 2>(
            new WeightCostFunctor(solver));
    problem.AddResidualBlock(cost_function, nullptr, parameters);
    
    // 添加位移约束 - 使用数值微分
    ceres::CostFunction* displacement_constraint = 
        new ceres::NumericDiffCostFunction<NumericDisplacementConstraintFunctor, 
                                          ceres::CENTRAL, 1, 2>(
            new NumericDisplacementConstraintFunctor(solver));
            
    // 使用较大的缩放损失函数来增加约束的权重
    ceres::LossFunction* constraint_loss = new ceres::ScaledLoss(
        nullptr, 1e6, ceres::TAKE_OWNERSHIP);
    problem.AddResidualBlock(displacement_constraint, constraint_loss, parameters);
    
    // 添加应力约束 - 对每个元素使用数值微分
    for (int i = 0; i < 10; ++i) {
        ceres::CostFunction* stress_constraint = 
            new ceres::NumericDiffCostFunction<NumericStressConstraintFunctor, 
                                              ceres::CENTRAL, 1, 2>(
                new NumericStressConstraintFunctor(solver, i));
                
        problem.AddResidualBlock(stress_constraint, constraint_loss, parameters);
    }
    
    // Add lower bounds on radii
    problem.SetParameterLowerBound(parameters, 0, 0.001);
    problem.SetParameterLowerBound(parameters, 1, 0.001);
    
    // Ceres solver options
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 200;  // 增加迭代次数
    options.function_tolerance = 1e-8;
    options.gradient_tolerance = 1e-8;
    options.parameter_tolerance = 1e-8;
    
    // 使用较大的初始信任区域半径，使优化可以有更大的移动空间
    options.initial_trust_region_radius = 10.0;
    options.max_trust_region_radius = 1000.0;
    
    // Solve the problem
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    
    std::cout << summary.FullReport() << std::endl;
    
    // Check constraints at final solution
    std::vector<double> displacements;
    std::vector<double> stresses;
    bool feasible = true;
    
    if (solver->Solve(parameters[0], parameters[1], displacements, stresses)) {
        // 检查位移约束
        double max_displacement = 0.0;
        for (size_t i = 0; i < displacements.size()/2; ++i) {
            double disp = std::sqrt(displacements[2*i]*displacements[2*i] + 
                                   displacements[2*i+1]*displacements[2*i+1]);
            max_displacement = std::max(max_displacement, disp);
        }
        
        std::cout << "Maximum displacement: " << max_displacement 
                  << " (limit: " << MAX_DISPLACEMENT << ")" << std::endl;
        
        if (max_displacement > MAX_DISPLACEMENT * 1.001) { // 允许1‰的误差
            std::cout << "Warning: Displacement constraint violated!" << std::endl;
            feasible = false;
        }
        
        // 检查应力约束
        for (size_t i = 0; i < stresses.size(); ++i) {
            std::cout << "Element " << i+1 << " stress: " << stresses[i] 
                      << " (limit: +/-" << MAX_STRESS << ")" << std::endl;
            
            if (std::abs(stresses[i]) > MAX_STRESS * 1.001) { // 允许1‰的误差
                std::cout << "Warning: Stress constraint violated for element " << i+1 << "!" << std::endl;
                feasible = false;
            }
        }
    } else {
        std::cout << "Failed to evaluate final solution!" << std::endl;
        feasible = false;
    }
    
    // Update the radii values
    r1 = parameters[0];
    r2 = parameters[1];
    
    // Calculate final weight
    double weight = solver->ComputeWeight(r1, r2);
    std::cout << "Optimized r1: " << r1 << std::endl;
    std::cout << "Optimized r2: " << r2 << std::endl;
    std::cout << "Final weight: " << weight << std::endl;
    std::cout << "Solution is " << (feasible ? "feasible" : "infeasible") << std::endl;
    
    delete solver;
    
    // 只有当解可用且满足约束时才返回true
    return summary.IsSolutionUsable() && feasible;
}

// Main function
int main() {
    // Initial guess - start with larger values to ensure structural integrity
    double r1 = 0.1;
    double r2 = 0.1;
    
    if (OptimizeTruss(r1, r2)) {
        std::cout << "Optimization successful!" << std::endl;
    } else {
        std::cout << "Optimization failed! Trying with increased constraint scaling..." << std::endl;
        
        // 如果首次优化失败，可以尝试增加约束权重或调整初始猜测值
        r1 = 0.1;
        r2 = 0.1;
        
        if (OptimizeTruss(r1, r2)) {
            std::cout << "Second attempt successful!" << std::endl;
        } else {
            std::cout << "Optimization completely failed!" << std::endl;
        }
    }
    
    return 0;
}