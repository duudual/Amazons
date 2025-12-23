#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cstring>

#define GRIDSIZE 8
#define INPUT_SIZE (GRIDSIZE * GRIDSIZE)  // 64个输入节点（棋盘状态）

using namespace std;

// 随机数生成器
random_device rd;
mt19937 gen(rd());
normal_distribution<double> normal_dist(0.0, 1.0);
uniform_real_distribution<double> uniform_dist(0.0, 1.0);

// 激活函数
class Activation {
public:
    virtual double forward(double x) = 0;
    virtual double backward(double x) = 0;
    virtual ~Activation() {}
};

class ReLU : public Activation {
public:
    double forward(double x) override {
        return max(0.0, x);
    }
    double backward(double x) override {
        return x > 0 ? 1.0 : 0.0;
    }
};

class Sigmoid : public Activation {
public:
    double forward(double x) override {
        return 1.0 / (1.0 + exp(-x));
    }
    double backward(double x) override {
        double s = forward(x);
        return s * (1.0 - s);
    }
};

class Tanh : public Activation {
public:
    double forward(double x) override {
        return tanh(x);
    }
    double backward(double x) override {
        double t = forward(x);
        return 1.0 - t * t;
    }
};

// 线性层
class Linear {
private:
    int input_size;
    int output_size;
    vector<vector<double>> weights;
    vector<double> biases;
    vector<vector<double>> weight_gradients;
    vector<double> bias_gradients;

public:
    Linear(int in_size, int out_size)
        : input_size(in_size), output_size(out_size) {
        // Xavier初始化权重
        double limit = sqrt(6.0 / (input_size + output_size));
        uniform_real_distribution<double> weight_dist(-limit, limit);

        weights.resize(output_size, vector<double>(input_size));
        biases.resize(output_size, 0.0);
        weight_gradients.resize(output_size, vector<double>(input_size, 0.0));
        bias_gradients.resize(output_size, 0.0);

        for (int i = 0; i < output_size; ++i) {
            for (int j = 0; j < input_size; ++j) {
                weights[i][j] = weight_dist(gen);
            }
        }
    }

    vector<double> forward(const vector<double>& input) {
        vector<double> output(output_size, 0.0);
        for (int i = 0; i < output_size; ++i) {
            for (int j = 0; j < input_size; ++j) {
                output[i] += weights[i][j] * input[j];
            }
            output[i] += biases[i];
        }
        return output;
    }

    vector<double> backward(const vector<double>& input, const vector<double>& output_grad) {
        // 计算输入梯度
        vector<double> input_grad(input_size, 0.0);
        for (int j = 0; j < input_size; ++j) {
            for (int i = 0; i < output_size; ++i) {
                input_grad[j] += weights[i][j] * output_grad[i];
            }
        }

        // 累积权重和偏置梯度
        for (int i = 0; i < output_size; ++i) {
            for (int j = 0; j < input_size; ++j) {
                weight_gradients[i][j] += output_grad[i] * input[j];
            }
            bias_gradients[i] += output_grad[i];
        }

        return input_grad;
    }

    void update_weights(double learning_rate, int batch_size) {
        for (int i = 0; i < output_size; ++i) {
            for (int j = 0; j < input_size; ++j) {
                weights[i][j] -= learning_rate * weight_gradients[i][j] / batch_size;
                weight_gradients[i][j] = 0.0;
            }
            biases[i] -= learning_rate * bias_gradients[i] / batch_size;
            bias_gradients[i] = 0.0;
        }
    }

    void save_weights(ofstream& file) {
        for (int i = 0; i < output_size; ++i) {
            for (int j = 0; j < input_size; ++j) {
                file << weights[i][j] << " ";
            }
            file << biases[i] << "\n";
        }
    }

    void load_weights(ifstream& file) {
        for (int i = 0; i < output_size; ++i) {
            for (int j = 0; j < input_size; ++j) {
                file >> weights[i][j];
            }
            file >> biases[i];
        }
    }
};

// 多层感知机价值网络
class ValueNetwork {
private:
    vector<Linear*> layers;
    vector<Activation*> activations;
    int num_layers;

public:
    ValueNetwork(const vector<int>& layer_sizes) {
        num_layers = layer_sizes.size() - 1;
        layers.resize(num_layers);
        activations.resize(num_layers - 1);  // 除了输出层外的所有层都有激活函数

        for (int i = 0; i < num_layers; ++i) {
            layers[i] = new Linear(layer_sizes[i], layer_sizes[i + 1]);

            if (i < num_layers - 1) {  // 隐藏层使用ReLU
                activations[i] = new ReLU();
            }
        }
    }

    ~ValueNetwork() {
        for (auto layer : layers) delete layer;
        for (auto act : activations) delete act;
    }

    // 将棋盘状态转换为网络输入
    vector<double> board_to_input(const int board[GRIDSIZE][GRIDSIZE]) {
        vector<double> input(INPUT_SIZE);
        for (int i = 0; i < GRIDSIZE; ++i) {
            for (int j = 0; j < GRIDSIZE; ++j) {
                int idx = i * GRIDSIZE + j;
                if (board[i][j] == 1) input[idx] = 1.0;      // 黑子
                else if (board[i][j] == -1) input[idx] = -1.0; // 白子
                else if (board[i][j] == 2) input[idx] = 0.5;   // 障碍
                else input[idx] = 0.0;                         // 空位
            }
        }
        return input;
    }

    // 前向传播
    double forward(const int board[GRIDSIZE][GRIDSIZE]) {
        vector<double> input = board_to_input(board);
        vector<double> current = input;

        for (int i = 0; i < num_layers; ++i) {
            current = layers[i]->forward(current);
            if (i < num_layers - 1) {  // 除了输出层外都要激活
                for (size_t j = 0; j < current.size(); ++j) {
                    current[j] = activations[i]->forward(current[j]);
                }
            }
        }

        // 输出层使用tanh激活，确保输出在[-1, 1]范围内
        return tanh(current[0]);
    }

    // 训练一步（单一样本）
    void train_step(const int board[GRIDSIZE][GRIDSIZE], double target, double learning_rate) {
        // 前向传播并保存中间结果
        vector<vector<double>> layer_inputs;
        vector<vector<double>> layer_outputs;

        vector<double> input = board_to_input(board);
        layer_inputs.push_back(input);
        vector<double> current = input;

        for (int i = 0; i < num_layers; ++i) {
            current = layers[i]->forward(current);
            layer_outputs.push_back(current);

            if (i < num_layers - 1) {
                for (size_t j = 0; j < current.size(); ++j) {
                    current[j] = activations[i]->forward(current[j]);
                }
            }
            layer_inputs.push_back(current);
        }

        // 计算输出层的损失梯度 (MSE损失对输出的导数)
        double output = layer_outputs.back()[0];
        double output_grad = 2.0 * (output - target) * (1.0 - output * output);  // tanh导数

        // 反向传播
        vector<double> current_grad = {output_grad};
        for (int i = num_layers - 1; i >= 0; --i) {
            if (i < num_layers - 1) {
                // 对激活函数求导
                for (size_t j = 0; j < current_grad.size(); ++j) {
                    current_grad[j] *= activations[i]->backward(layer_outputs[i][j]);
                }
            }
            current_grad = layers[i]->backward(layer_inputs[i], current_grad);
        }
    }

    // 更新权重（在batch结束后调用）
    void update_weights(double learning_rate, int batch_size) {
        for (auto layer : layers) {
            layer->update_weights(learning_rate, batch_size);
        }
    }

    // 保存模型
    void save_model(const string& filename) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Cannot open file " << filename << " for writing." << endl;
            return;
        }

        file << num_layers << "\n";
        for (auto layer : layers) {
            layer->save_weights(file);
        }
        file.close();
        cout << "Model saved to " << filename << endl;
    }

    // 加载模型
    void load_model(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Cannot open file " << filename << " for reading." << endl;
            return;
        }

        int saved_layers;
        file >> saved_layers;
        if (saved_layers != num_layers) {
            cerr << "Error: Model architecture mismatch." << endl;
            return;
        }

        for (auto layer : layers) {
            layer->load_weights(file);
        }
        file.close();
        cout << "Model loaded from " << filename << endl;
    }
};

// 训练数据结构
struct TrainingSample {
    int board[GRIDSIZE][GRIDSIZE];
    double target_value;

    TrainingSample() {
        memset(board, 0, sizeof(board));
        target_value = 0.0;
    }
};

// 批量训练函数
void train_network(ValueNetwork& network, const vector<TrainingSample>& samples,
                   int epochs, double learning_rate, int batch_size = 32) {
    cout << "Starting training with " << samples.size() << " samples..." << endl;

    for (int epoch = 0; epoch < epochs; ++epoch) {
        double total_loss = 0.0;
        vector<TrainingSample> batch;

        for (size_t i = 0; i < samples.size(); ++i) {
            batch.push_back(samples[i]);

            if (batch.size() == batch_size || i == samples.size() - 1) {
                // 训练一个batch
                for (const auto& sample : batch) {
                    network.train_step(sample.board, sample.target_value, learning_rate);
                }
                network.update_weights(learning_rate, batch.size());

                // 计算loss
                for (const auto& sample : batch) {
                    double pred = network.forward(sample.board);
                    total_loss += (pred - sample.target_value) * (pred - sample.target_value);
                }

                batch.clear();
            }
        }

        cout << "Epoch " << (epoch + 1) << "/" << epochs
             << ", Loss: " << total_loss / samples.size() << endl;
    }
}

// 测试函数
void test_network(ValueNetwork& network, const vector<TrainingSample>& test_samples) {
    double total_loss = 0.0;
    double total_abs_error = 0.0;

    for (const auto& sample : test_samples) {
        double pred = network.forward(sample.board);
        double error = pred - sample.target_value;
        total_loss += error * error;
        total_abs_error += abs(error);
    }

    cout << "Test Results:" << endl;
    cout << "MSE Loss: " << total_loss / test_samples.size() << endl;
    cout << "Mean Absolute Error: " << total_abs_error / test_samples.size() << endl;
}

int main() {
    // 创建一个3层网络: 64 -> 128 -> 64 -> 1
    vector<int> layer_sizes = {INPUT_SIZE, 128, 64, 1};
    ValueNetwork network(layer_sizes);

    // 示例：创建一些训练数据
    vector<TrainingSample> samples;

    // 这里可以加载实际的训练数据
    // 为了演示，我们创建一个简单的示例
    TrainingSample sample;
    // 设置一个简单的棋盘状态作为示例
    sample.board[0][0] = 1;  // 黑子
    sample.board[7][7] = -1; // 白子
    sample.target_value = 0.5; // 目标价值

    samples.push_back(sample);

    cout << "Network prediction before training: " << network.forward(sample.board) << endl;

    // 训练网络
    train_network(network, samples, 10, 0.01);

    cout << "Network prediction after training: " << network.forward(sample.board) << endl;

    // 保存模型
    network.save_model("value_network.model");

    // 加载模型测试
    ValueNetwork loaded_network(layer_sizes);
    loaded_network.load_model("value_network.model");

    cout << "Loaded network prediction: " << loaded_network.forward(sample.board) << endl;

    return 0;
}
