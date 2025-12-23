#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <ctime>

// 包含我们创建的头文件内容（为了演示，这里直接包含实现）
#include "value_network.cpp"
#include "data_generator.cpp"

int main() {
    srand(time(0));

    std::cout << "=== 亚马逊棋价值网络集成测试 ===\n\n";

    // 1. 生成一些训练数据（使用小参数以便快速测试）
    std::cout << "步骤1: 生成训练数据...\n";
    DataGenerator generator("test_mc_tree", 5, 5);  // max_depth=5, simulations=5
    generator.generate_training_data(1);  // 只生成1个游戏
    generator.save_training_data("test_training_data.bin");

    // 2. 加载训练数据
    std::cout << "\n步骤2: 加载训练数据...\n";
    auto training_data = generator.load_training_data("test_training_data.bin");

    if (training_data.empty()) {
        std::cout << "没有训练数据，使用合成数据...\n";
        // 创建一些合成训练数据
        TrainingSample sample1, sample2;
        // 简单局面1：只有黑子在角落
        sample1.board[0][0] = 1;
        sample1.mc_score = 0.3;

        // 简单局面2：黑白子都在角落
        sample2.board[0][0] = 1;
        sample2.board[7][7] = -1;
        sample2.mc_score = 0.0;

        training_data.push_back(sample1);
        training_data.push_back(sample2);
    }

    std::cout << "加载了 " << training_data.size() << " 个训练样本\n";

    // 3. 创建并训练价值网络
    std::cout << "\n步骤3: 创建和训练价值网络...\n";
    std::vector<int> layer_sizes = {INPUT_SIZE, 32, 16, 1};  // 64 -> 32 -> 16 -> 1
    ValueNetwork network(layer_sizes);

    // 测试训练前的预测
    std::cout << "训练前预测:\n";
    for (size_t i = 0; i < training_data.size(); ++i) {
        double pred = network.forward(training_data[i].board);
        std::cout << "样本 " << i << ": 真实值=" << training_data[i].mc_score
                  << ", 预测值=" << pred << "\n";
    }

    // 训练网络
    std::cout << "\n训练网络...\n";
    train_network(network, training_data, 20, 0.01, 1);  // 20个epoch，学习率0.01

    // 测试训练后的预测
    std::cout << "\n训练后预测:\n";
    double total_error = 0.0;
    for (size_t i = 0; i < training_data.size(); ++i) {
        double pred = network.forward(training_data[i].board);
        double error = abs(pred - training_data[i].mc_score);
        total_error += error;
        std::cout << "样本 " << i << ": 真实值=" << training_data[i].mc_score
                  << ", 预测值=" << pred << ", 误差=" << error << "\n";
    }

    std::cout << "\n平均绝对误差: " << total_error / training_data.size() << "\n";

    // 4. 保存和加载模型
    std::cout << "\n步骤4: 测试模型保存和加载...\n";
    network.save_model("test_value_network.model");

    ValueNetwork loaded_network(layer_sizes);
    loaded_network.load_model("test_value_network.model");

    // 验证加载的模型
    std::cout << "加载模型验证:\n";
    for (size_t i = 0; i < training_data.size(); ++i) {
        double original_pred = network.forward(training_data[i].board);
        double loaded_pred = loaded_network.forward(training_data[i].board);
        std::cout << "样本 " << i << ": 原模型=" << original_pred
                  << ", 加载模型=" << loaded_pred << "\n";
    }

    std::cout << "\n=== 测试完成 ===\n";
    std::cout << "生成的文件:\n";
    std::cout << "- test_training_data.bin (训练数据)\n";
    std::cout << "- test_value_network.model (训练好的模型)\n";
    std::cout << "- test_mc_tree_*.node (蒙特卡洛树节点文件)\n";

    return 0;
}



