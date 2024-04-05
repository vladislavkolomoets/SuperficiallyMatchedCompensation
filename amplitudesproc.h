#ifndef AMPLITUDESPROC_H
#define AMPLITUDESPROC_H


#include <unsupported/Eigen/FFT>

template<typename T>
void readFromBinFile(const std::string& _file_name, std::vector<std::vector<T>>& _matrix, std::pair<size_t, size_t>& _sizes_of_matrix)
{
    std::ifstream file(_file_name, std::ios::binary);

    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << _file_name << std::endl;
        return;
    }

    std::vector<T> current_vector{};
    T value{};

    size_t counter{};
    while (file.read(reinterpret_cast<char*>(&value), sizeof(T)))
    {
        current_vector.push_back(value);
        if (current_vector.size() == _sizes_of_matrix.second)
        {
            _matrix.push_back(current_vector);
            current_vector.clear();
        }
        counter++;

        if (counter > _sizes_of_matrix.first * _sizes_of_matrix.second)
        {
            break;
        }
    }

    file.close();
}

std::vector<std::vector<float>> transposeMatrix(const std::vector<std::vector<float>>& _matrix)
{
    std::vector<std::vector<float>> transposed_matrix{};

    if (_matrix.empty() || _matrix[0].empty())
    {
        std::cerr << "Vectors are empty." << std::endl;
        return transposed_matrix;
    }

    size_t num_rows = _matrix.size();
    size_t num_cols = _matrix[0].size();

    transposed_matrix.resize(num_cols, std::vector<float>(num_rows));

    for (size_t i = 0; i < num_rows; i++)
    {
        for (size_t j = 0; j < num_cols; j++)
        {
            transposed_matrix[j][i] = _matrix[i][j];
        }
    }

    return transposed_matrix;
}

void trimMatrix0(std::vector<std::vector<float>>& _matrix, size_t lower_bound, size_t upper_bound)
{
    if (lower_bound < 0 || upper_bound >= _matrix.size() || upper_bound < 0 || upper_bound >= _matrix.size() || lower_bound > upper_bound)
    {
        std::cerr << "Invalid bounds" << std::endl;
        return;
    }

    _matrix = transposeMatrix(_matrix);

    _matrix.erase(_matrix.begin(), _matrix.begin() + lower_bound);
    _matrix.erase(_matrix.begin() + (upper_bound - lower_bound + 1), _matrix.end());

    _matrix = transposeMatrix(_matrix);
}

void trimMatrix(std::vector<std::vector<float>>& _matrix, size_t lower_bound, size_t upper_bound)
{
    if (_matrix.size() == 0)
    {
        return;
    }

    if (lower_bound < 0 || upper_bound > _matrix[0].size())
    {
        return;
    }

    for (auto& row : _matrix)
    {
        row.erase(row.begin(), row.begin() + lower_bound);
        row.erase(row.begin() + (upper_bound - lower_bound + 1), row.end());
    }
}

std::vector<float> calculateRMS(const std::vector<std::vector<float>>& _matrix)
{
    std::vector<float> rms_values{};

    for (const auto& column : _matrix)
    {
        float sum_of_squares{};
        for (float value : column)
        {
            sum_of_squares += value * value;
        }
        float rms_current = std::sqrt(sum_of_squares / column.size());
        rms_values.push_back(rms_current);
    }

    return rms_values;
}

void writeVectorToFile(const std::vector<float>& _data, const std::string& _file_name)
{
    if (_data.size() == 0)
    {
        std::cerr << "The vector is empty\n";
    }

    //std::ofstream file(_file_name);
    std::ofstream file(_file_name, std::ios::out | std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << _file_name << std::endl;
        return;
    }

    /*
    for (float value : _data)
    {
        file << value << std::endl;
    }
    */

    file.write(reinterpret_cast<const char*>(_data.data()), _data.size() * sizeof(float));

    file.close();
}

Eigen::VectorXcf calculateFFT(const Eigen::VectorXf& _input)
{
    Eigen::FFT<float> fft{};

    Eigen::VectorXcf spectrum{};

    fft.fwd(spectrum, _input);

    return spectrum;
}

Eigen::VectorXf calculateAbs(const Eigen::VectorXcf& _spectrum)
{
    return _spectrum.cwiseAbs();
}

Eigen::VectorXf calculateIFFT(const Eigen::VectorXcf& _spectrum)
{
    Eigen::FFT<float> fft{};

    Eigen::VectorXf data{};

    fft.inv(data, _spectrum);

    return data;
}



std::vector<std::complex<float>> calculateFFT2(const std::vector<float>& _input)
{
    Eigen::FFT<float> fft{};

    Eigen::VectorXf input_eigen_vector = Eigen::Map<const Eigen::VectorXf>(_input.data(), _input.size());

    Eigen::VectorXcf spectrum{};

    fft.fwd(spectrum, input_eigen_vector);

    return std::vector<std::complex<float>>(spectrum.data(), spectrum.data() + spectrum.size());
}

std::vector<float> calculateAbs2(const std::vector<std::complex<float>>& _spectrum)
{
    Eigen::VectorXcf spectrum_eigen_vec(_spectrum.size());

    //for (size_t i = 0; i < spectrum_eigen_vec.size(); i++) // TODO : warning C4018: <: несоответствие типов со знаком и без знака
    for (Eigen::Index i = 0; i < spectrum_eigen_vec.size(); i++)
    {
        spectrum_eigen_vec(i) = _spectrum[i];
    }

    /*for (auto& spectrum_eigen_value : spectrum_eigen_vec)
    {
        spectrum_eigen_value = _spectrum[i];
    }*/

    Eigen::VectorXf amplitudes = spectrum_eigen_vec.cwiseAbs();

    return std::vector<float>(amplitudes.data(), amplitudes.data() + amplitudes.size());
}

std::vector<float> calculateIFFT2(const std::vector<std::complex<float>>& _spectrum)
{
    Eigen::FFT<float> fft{};

    Eigen::VectorXcf spectrum_eigen_vec(_spectrum.size());

    //for (size_t i = 0; i < spectrum_eigen_vec.size(); i++) // TODO : warning C4018: <: несоответствие типов со знаком и без знака
    for (Eigen::Index i = 0; i < spectrum_eigen_vec.size(); i++)
    {
        spectrum_eigen_vec(i) = _spectrum[i];
    }

    Eigen::VectorXf data{};

    fft.inv(data, spectrum_eigen_vec);

    return std::vector<float>(data.data(), data.data() + data.size());
}


Eigen::VectorXf vectorToEigen(const std::vector<float>& _vec)
{
    Eigen::VectorXf eigenVec(_vec.size());
    for (size_t i = 0; i < _vec.size(); i++)
    {
        eigenVec(i) = _vec[i];
    }
    return eigenVec;
}

std::vector<float> eigenToVector(const Eigen::VectorXf& _eigen_vec)
{
    std::vector<float> vec(_eigen_vec.size());
    for (int i = 0; i < _eigen_vec.size(); i++)
    {
        vec[i] = _eigen_vec(i);
    }
    return vec;
}

std::vector<std::vector<float>> f(std::vector<std::vector<float>>& _matrix)
{
    std::vector<std::vector<float>> result{};

    for (auto& vector : _matrix)
    {
        //auto temp_eigen_vec = vectorToEigen(vector);

        auto temp_fft = calculateFFT2(vector);

        auto temp_ampl = calculateAbs2(temp_fft);

        //auto temp_vec = eigenToVector(temp_ampl);

        result.push_back(temp_ampl);
    }

    return result;
}

Eigen::VectorXd convertFromVecFloatToEigenVecDouble(const std::vector<float>& _vector_float)
{
    Eigen::VectorXd vector_double(_vector_float.size());

    for (size_t i = 0; i < _vector_float.size(); ++i)
    {
        vector_double(i) = static_cast<double>(_vector_float[i]);
    }

    return vector_double;
}

#endif // LIBRARY_H