#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>
#include <random>
#define _USE_MATH_DEFINES
#include "math.h"
double function(double x)
{
    return sin(x) + 0.5;
}
double Evklidean_distance(double omega, double delta)
{
    return sqrt(omega * omega + delta * delta);
}
double Random(double min, double max)
{
    return (double)(rand()) / RAND_MAX * (max - min) + min;
}
std::vector<double> function_filtration(std::vector<double> f_noisy, std::vector <double> alpha)
{
    std::vector<double> result(f_noisy.size());
    int r = alpha.size();
    int M = r / 2;
    for (int k = 0; k < result.size(); k++)
    {
        result[k] = 0;
        for (int j = (k - M); j < (k + M + 1); j++)
        {
            if (j < 0 || j >(result.size() - 1))
                continue;
            result[k] += f_noisy[j] * alpha[j + M - k];
        }
    }
    return result;
}
double function_noisy(double x)
{
    double sigma = Random(-0.25, 0.25);
    return function(x) + sigma;
}
double omega(std::vector<double> f_filtation)//Критерий зашумленности
{
    double result = 0;
    double rez;
    for (int i = 1; i < f_filtation.size(); i++)
    {
        result += pow((f_filtation[i] - f_filtation[i - 1]), 2);
    }
    rez = sqrt(result);
    return rez;
}
double delta(std::vector<double> f_filtation, std::vector<double> f_noisy)//Критерий отличия
{
    double result = 0;
    double rez;
    for (int i = 0; i < f_filtation.size(); i++)
    {
        result += (pow((f_filtation[i] - f_noisy[i]), 2));
    }
    rez = sqrt((1. / f_filtation.size()) * result);
    return rez;
}
double J(double omega, double lamda, double delta)//
{
    return ((lamda * omega) + (1. - lamda) * delta);
}
std::vector<double> random_alpha(int r)
{
    double sum = 0;
    std::vector<double> result(r);
    for (int i = r / 2; i > 0; i--)
    {
        if (i == r / 2)
        {
            result[i] = Random(0, (1 - sum));
            sum += result[i];
        }
        else
        {
            result[i] = result[r - 1 - i] = Random(0, (1 - sum)) / 2;
            sum += result[i] + result[r - 1 - i];
        }
    }
    result[0] = result[r - 1] = (1 - sum) / 2;
    return result;
}
std::vector<double> random_search(std::vector<double> f_noisy, int r, double lamda)
{
    double e = 0.01;
    double P = 0.95;
    double Jmin = 0;
    int N = (log(1 - P) / log(1 - (e / M_PI)));
    std::vector<double> alpha;
    std::vector<double> f_filtration;
    std::vector<double> alpha_min;
    std::vector<double> result;
    double omega_min;
    double delta_min;
    double J_;
    double distance;
    for (int i = 0; i < N; i++)
    {
        alpha = random_alpha(r);
        f_filtration = function_filtration(f_noisy, alpha);
        J_ = J(omega(f_filtration), lamda, delta(f_filtration, f_noisy));
        if (i == 0 || Jmin > J_)
        {
            Jmin = J_;
            alpha_min = alpha;
            omega_min = omega(f_filtration);
            delta_min = delta(f_filtration, f_noisy);
        }
    }
    distance = Evklidean_distance(omega_min, delta_min);
    result.push_back(lamda);
    result.push_back(distance);
    for (int i = 0; i < r; i++)
    {
        result.push_back(alpha_min[i]);
    }
    result.push_back(omega_min);
    result.push_back(delta_min);
    result.push_back(Jmin);
    return result;
}
void print(std::vector<double> result, int r)
{
    for (int i = 0; i < result.size(); i++)
    {
        if (i == 0)
            std::cout << std::fixed << std::setprecision(1) << "|" << std::setw(3) << result[i] << "|";
        else
        {

            if (i > 1 && i < (2 + r))
            {
                if (i == 2)
                    std::cout << std::fixed << std::setprecision(4) << "[";
                if (i == 1 + r)
                    std::cout << std::fixed << std::setprecision(4) << std::setw(5) << result[i] << "]|";
                else
                    std::cout << std::fixed << std::setprecision(4) << std::setw(5) << result[i] << ", ";
            }
            else
                std::cout << std::fixed << std::setprecision(4) << std::setw(5) << result[i] << "|";
        }
    }
}
void print_(std::vector<double> result)
{
    int n = result.size();
    std::cout << "+-----+---------+--------+--------+" << std::endl;
    std::cout << "|  h* |    J    |    w   |   d    |" << std::endl;
    std::cout << "+-----+---------+--------+--------+" << std::endl;
    std::cout << std::fixed << std::setprecision(1) << "| " << result[0] << " |  " << std::fixed << std::setprecision(4) << result[n - 1] << " | " << result[n - 2] << " | " << result[n - 3] << " |" << std::endl;
    std::cout << "+-----+---------+--------+--------+" << std::endl;
}
void search(int r, std::vector<double> f_noisy)
{
    std::vector<double>result;
    std::vector<double>result_;
    double lamda;
    lamda = 0.0;
    double dist_min;
    if (r == 3)
    {
        std::cout << "+---+------+------------------------+------+------+------+" << std::endl;
        std::cout << "| h | dis  |          alpha         |  w   |  d   |  J   |" << std::endl;
        std::cout << "+---+------+------------------------+------+------+------+";
    }
    else
    {
        std::cout << "+---+------+----------------------------------------+------+------+------+" << std::endl;
        std::cout << "| h | dis  |                 alpha                  |  w   |  d   |  J   |" << std::endl;
        std::cout << "+---+------+----------------------------------------+------+------+------+";
    }
    for (int i = 0; i < 11; i++)
    {
        std::cout << std::endl;
        result = random_search(f_noisy, r, lamda);
        if (i == 0)
        {
            dist_min = result[1];
        }
        if (result[1] <= dist_min)
        {
            dist_min = result[1];
            result_ = result;
        }
        print(result, r);
        lamda = lamda + 0.1;
    }
    if (r == 3)
        std::cout << std::endl << "+---+------+------------------------+------+------+------+" << std::endl;
    else
        std::cout << std::endl << "+---+------+----------------------------------------+------+------+------+" << std::endl;
    print_(result_);

}
int main()
{
    srand(time(NULL));
    setlocale(LC_ALL, "rus");
    std::vector<double> f_noisy(101);
    for (int i = 0; i < 101; i++)
        f_noisy[i] = function_noisy((i * M_PI) / 100);
    std::cout << "Результаты численного эксперимента для r = 3 и оптимальное значение веса h*, функционала J и критериев w, d." << std::endl;
    search(3, f_noisy);
    std::cout << std::endl;
    std::cout << "\nРезультаты численного эксперимента для r = 5 и оптимальное значение веса h*, функционала J и критериев w, d." << std::endl;
    search(5, f_noisy);
    return 0;
}

