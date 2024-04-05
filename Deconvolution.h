#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H

#include <chrono>

#include "amplitudesproc.h"
#include "matrixproc.h"

class Deconvolution
{
private:

    std::vector<std::vector<float>> amplitudes_matrix{}; // untrimmed data with amplitudes

public:

    void processInputData(std::vector<std::string>&, std::pair<size_t, size_t>&, size_t, size_t, size_t, std::string);

    void calculateCorrections();

    void makeFiltration();

};

void Deconvolution::processInputData(
    std::vector<std::string>& _input_data_file_names,
    std::pair<size_t, size_t>& _sizes,
    size_t _dt,
    size_t _lower_bound,
    size_t _upper_bound,
    std::string _output_filename)
{
    /*=============================================================================================
                              I. Reading untrimmed data with amplitudes
    =============================================================================================*/

    auto start_reading_untrimmed_data = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tI. READING UNTRIMMED DATA WITH AMPLITUDES\n";

    std::string amplitudes_matrix_file_name = _input_data_file_names[0];

    readFromBinFile(amplitudes_matrix_file_name, this->amplitudes_matrix, _sizes); // float32

    std::cout << "Number of rows: " << amplitudes_matrix.size() << std::endl;

    auto end_reading_untrimmed_data = std::chrono::high_resolution_clock::now();

    auto duration_reading_untrimmed_data_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_reading_untrimmed_data - start_reading_untrimmed_data);

    std::cout << "\nReading untrimmed data time: \t" << duration_reading_untrimmed_data_ms.count()
        << " milliseconds" << std::endl;

    /*=============================================================================================
                                   II. Trim matrix with amplitudes
    =============================================================================================*/

    auto start_trim_matrix = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tII. TRIM MATRIX WITH AMPLITUDES\n";

    if (this->amplitudes_matrix.size() != 0)
    {
        trimMatrix(this->amplitudes_matrix, (_lower_bound / _dt), (_upper_bound / _dt));
    }

    std::cout << "After trim: " << this->amplitudes_matrix.size() << std::endl;

    auto end_trim_matrix = std::chrono::high_resolution_clock::now();

    auto duration_trim_matrix_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_trim_matrix - start_trim_matrix);

    std::cout << "\nTrim matrix time: \t" << duration_trim_matrix_ms.count()
        << " milliseconds" << std::endl;

    /*=============================================================================================
                           III. Calculate spectre of amplitudes
    =============================================================================================*/

    for (const auto& signal : this->amplitudes_matrix) 
    {
        std::vector<std::complex<float>> spectrum = calculateFFT2(signal);

        std::vector<float> amplitudes = calculateAbs2(spectrum);

        std::cout << "Amplitude Spectrum:" << std::endl;
        for (float amplitude : amplitudes) 
        {
            std::cout << amplitude << " ";
        }
        std::cout << std::endl;
    }
}

void Deconvolution::calculateCorrections()
{

}

void Deconvolution::makeFiltration()
{

}

#endif // DECONVOLUTION_H