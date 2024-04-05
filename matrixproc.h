#ifndef MATRIXPROC_H
#define MATRIXPROC_H

#include <iostream>
#include <array>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>
#include <memory>
#include <algorithm>
//
//#include <map>
#include <Eigen/Sparse>
//#include <Eigen/Dense>



using EnumerationNumber = std::map<int, size_t>;
using SpMatD = Eigen::SparseMatrix<double, Eigen::RowMajor>;
//using SpMatD = Eigen::SparseMatrix<double>;
//using SpMatUI = Eigen::SparseMatrix<size_t>;
//using SpMatF = Eigen::SparseMatrix<float>;

void writeToFile(const Eigen::VectorXd& _vector, const std::string _file_name)
{
    std::ofstream file(_file_name);
    if (!file.is_open())
    {
        std::cerr << "error with openning file: " << _file_name << std::endl;
        return;
    }

    for (auto& element : _vector)
    {
        file << element << "\n";
    }
    file.close();
}

void writeToFileWithNums(const Eigen::VectorXd& _vector, const std::vector<EnumerationNumber*>& _uniques_numbers, const std::vector<std::string>& _files_names, std::string _header = "")
{
    if (_uniques_numbers.size() != _files_names.size())
    {
        std::cerr << "_uniques_numbers.size() != _files_names.size()\n";
        return;
    }

    size_t shift = 0;

    for (size_t i = 0; i < _uniques_numbers.size(); ++i)
    {
        std::ofstream file(_header + _files_names[i], std::ios::out | std::ios::binary);
        if (!file.is_open())
        {
            std::cerr << "error with opening file: " << _files_names[i] << std::endl;
            return;
        }

        std::vector<std::pair<int, size_t>> sorted_uniques_numbers(_uniques_numbers[i]->begin(), _uniques_numbers[i]->end());

        std::sort(sorted_uniques_numbers.begin(), sorted_uniques_numbers.end(),
            [](const std::pair<int, size_t>& a, const std::pair<int, size_t>& b)
        {
            return a.second < b.second;
        });

        for (const auto& _unique : sorted_uniques_numbers)
        {
            file.write(reinterpret_cast<const char*>(&_unique.second), sizeof(_unique.second));
            file.write(reinterpret_cast<const char*>(&_vector[_unique.second + shift]), sizeof(double));
        }

        shift += _uniques_numbers[i]->size();

        file.close();
    }
}

Eigen::VectorXd readVectorFromBinaryFile(const std::string& _file_name, const size_t _rows_size)
{
    std::ifstream file(_file_name, std::ios::in | std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "error with openning file: " << _file_name << std::endl;
        return {};
    }

    double* data = new double[_rows_size];
    file.read((char*)data, sizeof(double) * _rows_size);

    Eigen::VectorXd result_Vector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(data, _rows_size);
    delete[] data;

    file.close();

    return result_Vector;
}

void logVec(Eigen::VectorXd& _vector)
{
    for (auto& element : _vector)
    {
        element = std::log(element);
    }
}

void substractAverage(Eigen::VectorXd& _vector)
{
    auto average = _vector.sum() / _vector.size();

    for (auto& element : _vector)
    {
        element -= average;
    }

}

void preproc(Eigen::VectorXd& _vector)
{
    auto average = _vector.sum() / _vector.size();

    int n_zeros = 0;

    for (int i = 0; i < _vector.size(); i++)
    {
        if (_vector[i] == 0)
        {
            n_zeros++;
        }
        _vector[i] = std::log(_vector[i]) - std::log(average);
    }

    std::cout << "Num_Zeros = " << n_zeros << "\n";
}

std::vector<Eigen::Triplet<size_t>> createTriplets(const std::vector<size_t>& _vector)
{
    if (_vector.size() == 0)
    {
        return{};
    }

    std::vector<Eigen::Triplet<size_t>> triplets{};

    for (size_t i = 0; i < _vector.size(); i++)
    {
        Eigen::Triplet<size_t> triplet(i, _vector[i], 1); // TODO : warning C4267: аргумент: преобразование из "size_t" в "const StorageIndex"; возможна потеря данных

        triplets.push_back(triplet);
    }

    return triplets;
}

void printTriplets(const std::vector<Eigen::Triplet<size_t>>& _triplets)
{
    //int i{};
    //int numbers_to_print = 20;

    for (const auto& triplet : _triplets)
    {
        //if (i < numbers_to_print)
        //{
        std::cout << "Triplet: (row = " << triplet.row() << ", col = " << triplet.col() << ", value = " << triplet.value() << ")" << std::endl;
        //}

        //i++;
    }
}

std::vector<int> findColsByRow(const std::vector<Eigen::Triplet<size_t>>& _triplets, size_t _n)
{
    std::vector<int> result{};

    //std::cout << "Row (" << _n << ") In cols: ";

    for (const auto& triplet : _triplets)
    {
        if (triplet.row() == _n)
        {
            //std::cout << triplet.col() << " ";
            result.push_back(triplet.col());
        }
    }

    //std::cout << "\n";
    return result;
}

Eigen::SparseMatrix<size_t> createMatrixFromTriplets(const std::vector<Eigen::Triplet<size_t>>& _triplets, const size_t _rows, const size_t _cols)
{
    if (_triplets.size() == 0 || _rows == 0 || _cols == 0)
    {
        return{};
    }

    Eigen::SparseMatrix<size_t> matrix(_rows, _cols);
    matrix.setFromTriplets(_triplets.begin(), _triplets.end());
    matrix.makeCompressed();
    return matrix;
}

std::vector<Eigen::Triplet<size_t>> concatenateTwoTriplets(const std::vector<Eigen::Triplet<size_t>>& _left_triplets, const std::vector<Eigen::Triplet<size_t>>& _right_triplets, const int _shift = 0)
{
    if (_left_triplets.size() == 0 && _right_triplets.size() != 0)
    {
        return _right_triplets;
    }
    else if (_left_triplets.size() != 0 && _right_triplets.size() == 0)
    {
        return _left_triplets;
    }
    else if (_left_triplets.size() == 0 && _right_triplets.size() == 0)
    {
        return {};
    }

    std::vector<Eigen::Triplet<size_t>> result_triplets(_left_triplets);

    for (size_t i = 0; i < _right_triplets.size(); i++)
    {
        Eigen::Triplet<size_t> temp_triplet = { _right_triplets.at(i).row(),
                                                _right_triplets.at(i).col() + _shift,
                                                _right_triplets.at(i).value() };

        result_triplets.push_back(temp_triplet);
        //std::cout << "Triplet: (x = " << temp_triplet.row() << ", y = " << temp_triplet.col() << ", value = " << temp_triplet.value() << ")" << std::endl;
    }

    return result_triplets;
}

std::vector<Eigen::Triplet<size_t>> concatenateTriplets(const std::vector<std::vector<Eigen::Triplet<size_t>>*>& _triplets, const std::vector<size_t>& _shifts)
{
    if (_triplets.size() != _shifts.size() || _triplets.size() == 0 || _shifts.size() == 0)
    {
        return {};
    }



    int shift = _shifts[0]; // TODO : warning C4267: инициализация: преобразование из "size_t" в "int"; возможна потеря данных
    auto result_triplet = concatenateTwoTriplets(*_triplets.at(0), *_triplets.at(1), shift);

    //std::vector<Eigen::Triplet<size_t>> result_triplet{};

    for (size_t i = 1; i < _triplets.size() - 1; i++)
    {
        //for (size_t j = 0; j < i; j++)
        //{
        shift += _shifts[i]; // TODO : warning C4267: +=: преобразование из "size_t" в "int"; возможна потеря данных
        //}

        result_triplet = concatenateTwoTriplets(result_triplet, *_triplets.at(i + 1), shift);
    }

    return result_triplet;
}

/*
    Calculate matrix equation: _A_matrix * result_Vector = _B_Vector with LSQR method
*/
Eigen::VectorXd solveLinearSystemWithLSQRMethod(const SpMatD& _A_matrix, const Eigen::VectorXd& _B_Vector)
{
    //Eigen::setNbThreads(8);
    //Eigen::initParallel();

    Eigen::LeastSquaresConjugateGradient<SpMatD> solver{};

    solver.setMaxIterations(200);

    solver.setTolerance(1e-5);

    //solver.compute(_A_matrix);

    auto start_analyze = std::chrono::high_resolution_clock::now();
    solver.analyzePattern(_A_matrix);
    auto end_analyze = std::chrono::high_resolution_clock::now();

    if (solver.info() != Eigen::Success)
    {
        std::cout << "LSQR failed\n";
        //auto duration_analyze_s = std::chrono::duration_cast<std::chrono::seconds>(end_analyze - start_analyze);
        //std::cout << "\n(in LSQR) Analyze time: \t" << duration_analyze_s.count() << " seconds" << std::endl;
        auto duration_analyze_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_analyze - start_analyze);
        std::cout << "(in LSQR) Analyze time: \t" << duration_analyze_ms.count() << " milliseconds" << std::endl;
        return {};
    }

    auto start_factorize = std::chrono::high_resolution_clock::now();
    solver.factorize(_A_matrix);
    auto end_factorize = std::chrono::high_resolution_clock::now();
    if (solver.info() != Eigen::Success)
    {
        std::cout << "LSQR failed\n";
        //auto duration_factorize_s = std::chrono::duration_cast<std::chrono::seconds>(end_factorize - start_factorize);
        //std::cout << "\n(in LSQR) Factorize time: \t" << duration_factorize_s.count() << " seconds" << std::endl;
        auto duration_factorize_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_factorize - start_factorize);
        std::cout << "(in LSQR) Factorize time: \t" << duration_factorize_ms.count() << " milliseconds" << std::endl;
        return {};
    }

    Eigen::VectorXd result_Vector = solver.solve(_B_Vector);

    std::cout << "#iterations:\t" << solver.iterations() << std::endl;
    std::cout << "estimated error:\t" << solver.error() << std::endl;

    //auto duration_analyze_s = std::chrono::duration_cast<std::chrono::seconds>(end_analyze - start_analyze);
    //std::cout << "\n(in LSQR) Analyze time: \t" << duration_analyze_s.count() << " seconds" << std::endl;
    auto duration_analyze_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_analyze - start_analyze);
    std::cout << "(in LSQR) Analyze time: \t" << duration_analyze_ms.count() << " milliseconds" << std::endl;

    //auto duration_factorize_s = std::chrono::duration_cast<std::chrono::seconds>(end_factorize - start_factorize);
    //std::cout << "\n(in LSQR) Factorize time: \t" << duration_factorize_s.count() << " seconds" << std::endl;
    auto duration_factorize_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_factorize - start_factorize);
    std::cout << "(in LSQR) Factorize time: \t" << duration_factorize_ms.count() << " milliseconds" << std::endl;

    return result_Vector;
    //return{};
}



/*
    Calculate matrix equation: _A_matrix * result_Vector = _B_Vector with Cholesky method
*/
Eigen::VectorXd solveLinearSystemWithCholesky(const SpMatD& _A_matrix, const Eigen::VectorXd& _B_Vector)
{
    //Eigen::setNbThreads(4);
    //Eigen::initParallel();

    // Cholesky decomposition
    Eigen::SimplicialLLT<SpMatD> solver{};

    //solver.compute(_A_matrix);

    auto start_analyze = std::chrono::high_resolution_clock::now();
    solver.analyzePattern(_A_matrix);
    auto end_analyze = std::chrono::high_resolution_clock::now();

    if (solver.info() != Eigen::Success)
    {
        std::cout << "Cholesky failed\n";
        //auto duration_analyze_s = std::chrono::duration_cast<std::chrono::seconds>(end_analyze - start_analyze);
        //std::cout << "\n(in LSQR) Analyze time: \t" << duration_analyze_s.count() << " seconds" << std::endl;
        auto duration_analyze_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_analyze - start_analyze);
        std::cout << "(in Cholesky) Analyze time: \t" << duration_analyze_ms.count() << " milliseconds" << std::endl;
        return {};
    }

    auto start_factorize = std::chrono::high_resolution_clock::now();
    solver.factorize(_A_matrix);
    auto end_factorize = std::chrono::high_resolution_clock::now();
    if (solver.info() != Eigen::Success)
    {
        std::cout << "Cholesky failed\n";
        //auto duration_factorize_s = std::chrono::duration_cast<std::chrono::seconds>(end_factorize - start_factorize);
        //std::cout << "\n(in LSQR) Factorize time: \t" << duration_factorize_s.count() << " seconds" << std::endl;
        auto duration_factorize_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_factorize - start_factorize);
        std::cout << "(in Cholesky) Factorize time: \t" << duration_factorize_ms.count() << " milliseconds" << std::endl;
        return {};
    }

    std::cout << "Cholesky decomposition succeeded\n";

    // Solve using the Cholesky decomposition
    Eigen::VectorXd result_Vector = solver.solve(_B_Vector);

    if (solver.info() != Eigen::Success)
    {
        std::cout << "Cholesky solver failed\n";
        return {};
    }

    std::cout << "Cholesky solver succeeded\n";

    //auto duration_analyze_s = std::chrono::duration_cast<std::chrono::seconds>(end_analyze - start_analyze);
    //std::cout << "\n(in LSQR) Analyze time: \t" << duration_analyze_s.count() << " seconds" << std::endl;
    auto duration_analyze_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_analyze - start_analyze);
    std::cout << "(in Cholesky) Analyze time: \t" << duration_analyze_ms.count() << " milliseconds" << std::endl;

    //auto duration_factorize_s = std::chrono::duration_cast<std::chrono::seconds>(end_factorize - start_factorize);
    //std::cout << "\n(in LSQR) Factorize time: \t" << duration_factorize_s.count() << " seconds" << std::endl;
    auto duration_factorize_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_factorize - start_factorize);
    std::cout << "(in Cholesky) Factorize time: \t" << duration_factorize_ms.count() << " milliseconds" << std::endl;


    return result_Vector;
}


// https://eigen.tuxfamily.org/dox/classEigen_1_1BiCGSTAB.html
Eigen::VectorXd solveLinearSystemWithBICGSTAB(const SpMatD& _A_matrix, const Eigen::VectorXd& _B_Vector)
{
    //Eigen::setNbThreads(4);
    //Eigen::initParallel();

    Eigen::BiCGSTAB<SpMatD> solver{};

    solver.setTolerance(1e-6);

    solver.setMaxIterations(500);

    //solver.compute(_A_matrix);

    auto start_analyze = std::chrono::high_resolution_clock::now();
    solver.analyzePattern(_A_matrix);
    auto end_analyze = std::chrono::high_resolution_clock::now();

    if (solver.info() != Eigen::Success)
    {
        std::cout << "BICGSTAB failed\n";
        //auto duration_analyze_s = std::chrono::duration_cast<std::chrono::seconds>(end_analyze - start_analyze);
        //std::cout << "\n(in LSQR) Analyze time: \t" << duration_analyze_s.count() << " seconds" << std::endl;
        auto duration_analyze_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_analyze - start_analyze);
        std::cout << "(in BICGSTAB) Analyze time: \t" << duration_analyze_ms.count() << " milliseconds" << std::endl;
        return {};
    }

    auto start_factorize = std::chrono::high_resolution_clock::now();
    solver.factorize(_A_matrix);
    auto end_factorize = std::chrono::high_resolution_clock::now();
    if (solver.info() != Eigen::Success)
    {
        std::cout << "BICGSTAB failed\n";
        //auto duration_factorize_s = std::chrono::duration_cast<std::chrono::seconds>(end_factorize - start_factorize);
        //std::cout << "\n(in LSQR) Factorize time: \t" << duration_factorize_s.count() << " seconds" << std::endl;
        auto duration_factorize_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_factorize - start_factorize);
        std::cout << "(in BICGSTAB) Factorize time: \t" << duration_factorize_ms.count() << " milliseconds" << std::endl;
        return {};
    }

    Eigen::VectorXd result_Vector = solver.solve(_B_Vector);

    //auto duration_analyze_s = std::chrono::duration_cast<std::chrono::seconds>(end_analyze - start_analyze);
    //std::cout << "\n(in LSQR) Analyze time: \t" << duration_analyze_s.count() << " seconds" << std::endl;
    auto duration_analyze_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_analyze - start_analyze);
    std::cout << "(in BICGSTAB) Analyze time: \t" << duration_analyze_ms.count() << " milliseconds" << std::endl;

    //auto duration_factorize_s = std::chrono::duration_cast<std::chrono::seconds>(end_factorize - start_factorize);
    //std::cout << "\n(in LSQR) Factorize time: \t" << duration_factorize_s.count() << " seconds" << std::endl;
    auto duration_factorize_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_factorize - start_factorize);
    std::cout << "(in BICGSTAB) Factorize time: \t" << duration_factorize_ms.count() << " milliseconds" << std::endl;

    std::cout << "#iterations:\t" << solver.iterations() << std::endl;
    std::cout << "estimated error:\t" << solver.error() << std::endl;

    return result_Vector;
}

#endif // MATRIXPROC_H