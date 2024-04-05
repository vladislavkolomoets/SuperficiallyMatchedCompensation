#ifndef AMPLITUDES_CORRECTION_H
#define AMPLITUDES_CORRECTION_H

#include <chrono>

#include "dataproc.h"
#include "matrixproc.h"
#include "amplitudesproc.h"
#include "correctproc.h"

const int ROWS_SIZE = 6147792; //71969600 - 1; //12;

template<typename T>
class SuperficiallyMatchedCompensation
{
private:

    std::vector<std::vector<T>> amplitudes_matrix{}; // untrimmed data with amplitudes

    Eigen::VectorXd amplitudes{}; // amplitudes after rms

    /* numbers of factors */
    std::unique_ptr<int[]> sou_numbers{};
    std::unique_ptr<int[]> rec_numbers{};
    std::unique_ptr<int[]> cdp_numbers{};
    std::unique_ptr<int[]> off_numbers{};
    std::unique_ptr<int[]> az_numbers{};

    size_t n_factors{}; // actual number of factors

    /* dictionary of unique numbers of factors */
    std::unique_ptr<EnumerationNumber> sou_uniques{};
    std::unique_ptr<EnumerationNumber> rec_uniques{};
    std::unique_ptr<EnumerationNumber> cdp_uniques{};
    std::unique_ptr<EnumerationNumber> off_uniques{};
    std::unique_ptr<EnumerationNumber> az_uniques{};

    /* input numbers after coding by dictionary (F_uniques) */
    std::unique_ptr<std::vector<size_t>> sou_codes{};
    std::unique_ptr<std::vector<size_t>> rec_codes{};
    std::unique_ptr<std::vector<size_t>> cdp_codes{};
    std::unique_ptr<std::vector<size_t>> off_codes{};
    std::unique_ptr<std::vector<size_t>> az_codes{};

    /* triplets containing positions of '1'-s for each of factors */
    std::unique_ptr<std::vector<Eigen::Triplet<size_t>>> sou_triplets{};
    std::unique_ptr<std::vector<Eigen::Triplet<size_t>>> rec_triplets{};
    std::unique_ptr<std::vector<Eigen::Triplet<size_t>>> cdp_triplets{};
    std::unique_ptr<std::vector<Eigen::Triplet<size_t>>> off_triplets{};
    std::unique_ptr<std::vector<Eigen::Triplet<size_t>>> az_triplets{};

    std::vector<Eigen::Triplet<size_t>> concatenated_triplets{}; // triplets after merging of all locals triplets

    std::unique_ptr<Eigen::SparseMatrix<size_t>> conc_mat{}; // matrix created from concatenated_triplets

    int obs_mat_sys_rows{}; // number of rows in observe system matrix

    Eigen::VectorXd solution{}; // solution of equation A x = b

    double amplitude_average{}; // average value of amplitudes vector

    std::vector<double> corrections{}; // amplitudes processing with corrections

public:

    void processInputData(std::vector<std::string>&, std::pair<size_t, size_t>&, size_t, size_t, size_t, std::string);

    void calculateCorrections();

    void enterCorrections();
};

template<typename T>
void SuperficiallyMatchedCompensation<T>::processInputData(
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

    readFromBinFile<T>(amplitudes_matrix_file_name, this->amplitudes_matrix, _sizes); // float32

    std::cout << "Number of rows: " << amplitudes_matrix.size() << std::endl;

    std::cout << "+++++++\n";
    for (size_t i = 0; i < 5; i++)
    {
        for (size_t j = 0; j < 100; j++)
        {
            std::cout << this->amplitudes_matrix[i][j] << " ";
        }

        std::cout << std::endl;
    }
    std::cout << "+++++++\n";

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

    trimMatrix(this->amplitudes_matrix, (_lower_bound / _dt), (_upper_bound / _dt));

    std::cout << "After trim: " << this->amplitudes_matrix.size() << std::endl;

    std::cout << "+++++++\n";
    for (size_t i = 0; i < 5; i++)
    {
        for (size_t j = 0; j < 100; j++)
        {
            std::cout << this->amplitudes_matrix[i][j] << " ";
        }

        std::cout << std::endl;
    }
    std::cout << "+++++++\n";

    auto end_trim_matrix = std::chrono::high_resolution_clock::now();

    auto duration_trim_matrix_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_trim_matrix - start_trim_matrix);

    std::cout << "\nTrim matrix time: \t" << duration_trim_matrix_ms.count()
        << " milliseconds" << std::endl;

    /*=============================================================================================
                           III. Calculate RMS for trimmed matrix of amplitudes
    =============================================================================================*/

    auto start_calc_rms = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tIII. CALCULATE RMS FOR TRIMMED MATRIX OF AMPLITUDES\n";

    std::vector<float> rms_values = calculateRMS(this->amplitudes_matrix);

    std::cout << "RMS values: " << rms_values.size() << std::endl;

    std::cout << "square\n";
    for (size_t j = 0; j < 100; j++)
    {
        std::cout << this->amplitudes_matrix[4669400][j] << " ";
    }
    std::cout << "\n"; int sum{};
    for (size_t j = 0; j < 100; j++)
    {
        std::cout << (this->amplitudes_matrix[4669400][j]) * (this->amplitudes_matrix[4669400][j]) << " ";
        sum += (this->amplitudes_matrix[4669400][j]) * (this->amplitudes_matrix[4669400][j]);
    }
    std::cout << "\nsum = " << sum << "\n";
    std::cout << "\nsize = " << amplitudes_matrix[4669400].size() << "\n";
    std::cout << "\nsqrt = " << std::sqrt(sum / amplitudes_matrix[4669400].size()) << "\n";
    for (size_t j = 0; j < this->amplitudes_matrix[4669400].size(); j++)
    {
        std::cout << this->amplitudes_matrix[4669400][j] << " ";
    }
    std::cout << std::endl;

    this->amplitudes = convertFromVecFloatToEigenVecDouble(rms_values);

    this->amplitudes_matrix.clear();

    auto end_calc_rms = std::chrono::high_resolution_clock::now();

    auto duration_calc_rms_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_calc_rms - start_calc_rms);

    std::cout << "\nCalculate RMS time: \t" << duration_calc_rms_ms.count()
        << " milliseconds" << std::endl;

    auto start_write_rms_to_file = std::chrono::high_resolution_clock::now();

    writeVectorToFile(rms_values, _output_filename);

    //rms_values.clear();

    auto end_write_rms_to_file = std::chrono::high_resolution_clock::now();

    auto duration_write_rms_to_file_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_write_rms_to_file - start_write_rms_to_file);

    std::cout << "\nWrite RMS ro file time: \t" << duration_write_rms_to_file_ms.count()
        << " milliseconds" << std::endl;

    /*=============================================================================================
                           IV. Read input files with numbers
    =============================================================================================*/

    auto start_reading_data = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tIV. CALCULATE RMS FOR TRIMMED MATRIX OF AMPLITUDES\n";

    std::string sou_file_name = _input_data_file_names[1];
    std::string rec_file_name = _input_data_file_names[2];
    std::string cdp_file_name = _input_data_file_names[3];
    std::string off_file_name = _input_data_file_names[4];
    std::string az_file_name = _input_data_file_names[5];

    /*int* sou_numbers = nullptr;
    int* rec_numbers = nullptr;
    int* cdp_numbers = nullptr;
    int* off_numbers = nullptr;
    int* az_numbers = nullptr;*/


    this->sou_numbers = readColomnFromBinaryFile(sou_file_name, ROWS_SIZE);
    if (this->sou_numbers != nullptr)
    {
        this->n_factors++;
    }

    this->rec_numbers = readColomnFromBinaryFile(rec_file_name, ROWS_SIZE);
    if (this->rec_numbers != nullptr)
    {
        this->n_factors++;
    }

    this->cdp_numbers = readColomnFromBinaryFile(cdp_file_name, ROWS_SIZE);
    if (this->cdp_numbers != nullptr)
    {
        this->n_factors++;
    }

    this->off_numbers = readColomnFromBinaryFile(off_file_name, ROWS_SIZE);
    if (this->off_numbers != nullptr)
    {
        this->n_factors++;
    }

    this->az_numbers = readColomnFromBinaryFile(az_file_name, ROWS_SIZE);
    if (this->az_numbers != nullptr)
    {
        this->n_factors++;
    }

    auto end_reading_data = std::chrono::high_resolution_clock::now();

    auto duration_reading_data_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_reading_data - start_reading_data);

    std::cout << "\nReading data time: \t" << duration_reading_data_ms.count()
        << " milliseconds" << std::endl;

    /*=============================================================================================
                                        V. Enumerate
    =============================================================================================*/

    auto start_calc_enumeration = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tV. ENUMERATING\n";

    this->sou_uniques = std::make_unique<EnumerationNumber>(calcEnumerationNumber(this->sou_numbers.get(), ROWS_SIZE));
    this->rec_uniques = std::make_unique<EnumerationNumber>(calcEnumerationNumber(this->rec_numbers.get(), ROWS_SIZE));
    this->cdp_uniques = std::make_unique<EnumerationNumber>(calcEnumerationNumber(this->cdp_numbers.get(), ROWS_SIZE));
    this->off_uniques = std::make_unique<EnumerationNumber>(calcEnumerationNumber(this->off_numbers.get(), ROWS_SIZE));
    this->az_uniques = std::make_unique<EnumerationNumber>(calcEnumerationNumber(this->az_numbers.get(), ROWS_SIZE));

    auto end_calc_enumeration = std::chrono::high_resolution_clock::now();

    auto duration_calc_enumeration_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_calc_enumeration - start_calc_enumeration);

    std::cout << "\nCalc enums time: \t" << duration_calc_enumeration_ms.count()
        << " milliseconds" << std::endl;

    /*=============================================================================================
                                        VI. Coding
    =============================================================================================*/

    auto start_make_codes = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tVI. CODING\n";

    this->sou_codes = std::make_unique<std::vector<size_t>>(makeCodes(*(this->sou_uniques), this->sou_numbers.get(), ROWS_SIZE));
    this->rec_codes = std::make_unique<std::vector<size_t>>(makeCodes(*(this->rec_uniques), this->rec_numbers.get(), ROWS_SIZE));
    this->cdp_codes = std::make_unique<std::vector<size_t>>(makeCodes(*(this->cdp_uniques), this->cdp_numbers.get(), ROWS_SIZE));
    this->off_codes = std::make_unique<std::vector<size_t>>(makeCodes(*(this->off_uniques), this->off_numbers.get(), ROWS_SIZE));
    this->az_codes = std::make_unique<std::vector<size_t>>(makeCodes(*(this->az_uniques), this->az_numbers.get(), ROWS_SIZE));

    /*auto* sou_codes =
        new std::vector<size_t>(makeCodes(*sou_uniques, sou_numbers.get(), ROWS_SIZE));
    auto* rec_codes =
        new std::vector<size_t>(makeCodes(*rec_uniques, rec_numbers.get(), ROWS_SIZE));
    auto* cdp_codes =
        new std::vector<size_t>(makeCodes(*cdp_uniques, cdp_numbers.get(), ROWS_SIZE));
    auto* off_codes =
        new std::vector<size_t>(makeCodes(*off_uniques, off_numbers.get(), ROWS_SIZE));
    auto* az_codes =
        new std::vector<size_t>(makeCodes(*az_uniques, az_numbers.get(), ROWS_SIZE));*/

        //delete[] sou_numbers;
    this->sou_numbers.reset();
    //delete[] rec_numbers;
    this->rec_numbers.reset();
    //delete[] cdp_numbers;
    this->cdp_numbers.reset();
    //delete[] off_numbers;
    this->off_numbers.reset();
    //delete[] az_numbers;
    this->az_numbers.reset();

    auto end_make_codes = std::chrono::high_resolution_clock::now();

    auto duration_make_codes_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_make_codes - start_make_codes);

    std::cout << "\nMaking codes time: \t" << duration_make_codes_ms.count()
        << " milliseconds" << std::endl;
}

template<typename T>
void SuperficiallyMatchedCompensation<T>::calculateCorrections()
{
    /*=============================================================================================
                                        I. Creating triplets
    =============================================================================================*/

    auto start_creating_triplets = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tI. CREATING TRIPLETS\n";

    this->sou_triplets = std::make_unique<std::vector<Eigen::Triplet<size_t>>>(createTriplets(*(this->sou_codes)));
    this->rec_triplets = std::make_unique<std::vector<Eigen::Triplet<size_t>>>(createTriplets(*(this->rec_codes)));
    this->cdp_triplets = std::make_unique<std::vector<Eigen::Triplet<size_t>>>(createTriplets(*(this->cdp_codes)));
    this->off_triplets = std::make_unique<std::vector<Eigen::Triplet<size_t>>>(createTriplets(*(this->off_codes)));
    this->az_triplets = std::make_unique<std::vector<Eigen::Triplet<size_t>>>(createTriplets(*(this->az_codes)));

    /*auto* sou_triplets =
        new std::vector<Eigen::Triplet<size_t>>(createTriplets(*sou_codes));
    auto* rec_triplets =
        new std::vector<Eigen::Triplet<size_t>>(createTriplets(*rec_codes));
    auto* cdp_triplets =
        new std::vector<Eigen::Triplet<size_t>>(createTriplets(*cdp_codes));
    auto* off_triplets =
        new std::vector<Eigen::Triplet<size_t>>(createTriplets(*off_codes));
    auto* az_triplets =
        new std::vector<Eigen::Triplet<size_t>>(createTriplets(*az_codes));*/

    std::vector<std::vector<Eigen::Triplet<size_t>>*> triplets =
    {
        this->sou_triplets.get(),
        this->rec_triplets.get(),
        this->cdp_triplets.get(),
        this->off_triplets.get(),
        this->az_triplets.get()
    };

    std::vector<size_t> shifts =
    {
        this->sou_uniques->size(),
        this->rec_uniques->size(),
        this->cdp_uniques->size(),
        this->off_uniques->size(),
        this->az_uniques->size()
    };


    this->concatenated_triplets = concatenateTriplets(triplets, shifts);

    std::cout << "\nTriplets size = " << concatenated_triplets.size() << std::endl;

    /*std::cout << "sou_triplets:\n";
    printTriplets(*sou_triplets);
    std::cout << "rec_triplets:\n";
    printTriplets(*rec_triplets);
    std::cout << "concatenated_triplets:\n";
    printTriplets(concatenated_triplets);*/

    // std::cout << SOU_triplets->size() << "\n" << REC_triplets->size() << "\n"
    // << concatenated_triplets.size() << "\n"; printTriplets(*SOU_triplets);
    // std::cout << "\n"; printTriplets(*REC_triplets); std::cout << "\n";
    // printTriplets(concatenated_triplets); std::cout << "\n";

    //delete sou_triplets;
    this->sou_triplets.reset();
    //delete rec_triplets;
    this->rec_triplets.reset();
    //delete cdp_triplets;
    this->cdp_triplets.reset();
    //delete off_triplets;
    this->off_triplets.reset();
    //delete az_triplets;
    this->az_triplets.reset();
    //delete sou_codes;
    this->sou_codes.reset();
    //delete rec_codes;
    this->rec_codes.reset();
    //delete cdp_codes;
    this->cdp_codes.reset();
    //delete off_codes;
    this->off_codes.reset();
    //delete az_codes;
    this->az_codes.reset();

    auto end_creating_triplets = std::chrono::high_resolution_clock::now();

    auto duration_creating_triplets_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_creating_triplets - start_creating_triplets);

    std::cout << "\nCreating tr-s time: \t"
        << duration_creating_triplets_ms.count() << " milliseconds" << std::endl;

    /*=============================================================================================
                                        II. Creating Mat
    =============================================================================================*/

    auto start_creating_mats = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tII. CREATING MATS\n";

    this->conc_mat = std::make_unique<Eigen::SparseMatrix<size_t>>(createMatrixFromTriplets
    (
        this->concatenated_triplets, ROWS_SIZE,
        this->sou_uniques->size() +
        this->rec_uniques->size() +
        this->cdp_uniques->size() +
        this->off_uniques->size() +
        this->az_uniques->size()
    ));

    /*auto* conc_mat = new Eigen::SparseMatrix<size_t>(createMatrixFromTriplets
    (
        concatenated_triplets, ROWS_SIZE,
        sou_uniques->size() +
        rec_uniques->size() +
        cdp_uniques->size() +
        off_uniques->size() +
        az_uniques->size()
    ));*/


    auto end_creating_mats = std::chrono::high_resolution_clock::now();

    auto duration_creating_mats_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_creating_mats - start_creating_mats);

    std::cout << "\nCreating mats time: \t" << duration_creating_mats_ms.count()
        << " milliseconds" << std::endl;

    /*=============================================================================================
                                        III. Reading amplitudes
    =============================================================================================*/

    auto start_reading_ampls = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tIII. READING AMPLITUDES\n";

    //logVec(amplitudes); substractAverage(amplitudes);

    this->amplitude_average = this->amplitudes.sum() / this->amplitudes.size();

    preproc(this->amplitudes);

    std::cout << "+++++++\n";
    for (size_t i = 0; i < 20; i++)
    {
        std::cout << this->amplitudes(i) << " ";

        std::cout << std::endl;
    }
    std::cout << "+++++++\n";

    auto end_reading_ampls = std::chrono::high_resolution_clock::now();

    auto duration_reading_ampls_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_reading_ampls - start_reading_ampls);

    std::cout << "\nReading ampls time: \t" << duration_reading_ampls_ms.count()
        << " milliseconds" << std::endl;

    /*=============================================================================================
                                        IV. Solving system
    =============================================================================================*/

    auto start_solving_sys = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tIV. SOLVING SYS\n";

    SpMatD obs_sys_mat = this->conc_mat->cast<double>();
    //delete conc_mat;
    this->conc_mat.reset();

    /*LSQR*/
    // this->solution = solveLinearSystemWithLSQRMethod(obs_sys_mat, amplitudes);

    /*LU*/
    // auto obs_sys_mat_transpose = obs_sys_mat.transpose();
    // SpMatF B = obs_sys_mat_transpose * obs_sys_mat;
    // Eigen::VectorXf y = obs_sys_mat_transpose * (amplitudes);
    // this->solution = solveLinearSystemWithLU(B, y);

    /*Cholesky*/
    // auto obs_sys_mat_transpose = obs_sys_mat.transpose();
    // SpMatD B = obs_sys_mat_transpose * obs_sys_mat;
    // Eigen::VectorXd y = obs_sys_mat_transpose * (amplitudes);
    // this->solution = solveLinearSystemWithCholesky(obs_sys_mat, amplitudes);
    // obs_sys_mat.resize(0, 0);
    // obs_sys_mat_transpose.resize(0, 0);

    /*QR*/
    // this->solution = solveLinearSystemWithQR(obs_sys_mat, amplitudes);

    /*BICGSTAB*/
    //auto start_rect_to_square = std::chrono::high_resolution_clock::now();
    auto obs_sys_mat_transpose = obs_sys_mat.transpose();
    SpMatD B = obs_sys_mat_transpose * obs_sys_mat;
    Eigen::VectorXd y = obs_sys_mat_transpose * (amplitudes); writeToFile(y, "ddd");
    // auto end_rect_to_square = std::chrono::high_resolution_clock::now();

    this->solution = solveLinearSystemWithBICGSTAB(B, y);

    this->obs_mat_sys_rows = obs_sys_mat.rows();

    std::cout << "Matrix.size = " << obs_sys_mat.size()
        << "\tMatrix.rows = " << obs_sys_mat.rows()
        << "\tMatrix.cols = " << obs_sys_mat.cols() << "\n";
    std::cout << "Amplitudes.size = " << amplitudes.size() << "\n";
    std::cout << "Solution.size = " << solution.size() << "\n";

    //auto X = obs_sys_mat * solution;

    //std::cout << X.rows() << "\t" << X.cols() << "\n" << X << "\n";

    obs_sys_mat.resize(0, 0);
    obs_sys_mat_transpose.resize(0, 0);

    auto end_solving_sys = std::chrono::high_resolution_clock::now();

    auto duration_solving_sys_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_solving_sys - start_solving_sys);

    std::cout << "\nSolving sys time: \t" << duration_solving_sys_ms.count()
        << " milliseconds" << std::endl;

    /*=============================================================================================
                                       V. Writing to file
    =============================================================================================*/

    auto start_writing_out = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tV. WRITING OUT\n";

    // writeToFile(solution,
    // "/mnt/c/users/geouser7/source/repos/ConsoleApplication1/ConsoleApplication1/files2/lsqr_solution.txt");

    std::vector<EnumerationNumber*> uniques_numbers =
    {
        this->sou_uniques.get(),
        this->rec_uniques.get(),
        this->cdp_uniques.get(),
        this->off_uniques.get(),
        this->az_uniques.get()
    };

    std::vector<std::string> files_names =
    {
        "sou.bin",
        "rec.bin",
        "cdp.bin",
        "off.bin",
        "az.bin"
    };

    std::string header = "solution_lsqr_"; // lsqr // BICGSTAB // Cholesky
    //std::cout << "\nSolution:\n" << solution << "\n";
    //std::cout << "\nAmplitudes:\n" << amplitudes << "\n";
    //std::cout << "\*****:\n" << amplitudes * solution.transpose() << "\n";
    writeToFileWithNums(this->solution, uniques_numbers, files_names, header);


    auto end_writing_out = std::chrono::high_resolution_clock::now();

    auto duration_writing_out_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_writing_out - start_writing_out);

    std::cout << "\nWriting out time: \t" << duration_writing_out_ms.count()
        << " milliseconds" << std::endl;
}

template<typename T>
void SuperficiallyMatchedCompensation<T>::enterCorrections()
{ // TODO : message : во врем€ компил€ции функции-члена класс шаблон "void SuperficiallyMatchedCompensation<float>::enterCorrections(void)"
    /*=============================================================================================
                            I. Entering corrections into amplitude vector
    =============================================================================================*/

    auto start_entering_corrections = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tI. ENTERING CORRECTIONS INTO AMPLITUDE VECTOR\n";

    //std::cout << "matrix:\n" << obs_sys_mat << std::endl;
    //std::cout << "solution:\n" << solution << std::endl;
    //std::cout << "amplitudes:\n" << amplitudes << std::endl;

    std::cout << "obs_mat_sys_rows:\n" << obs_mat_sys_rows << std::endl;

    /*std::vector<std::vector<int>> vector_of_global_cols{};
    for (size_t i = 0; i < this->obs_mat_sys_rows; i++)
    {
        auto v = findColsByRow(this->concatenated_triplets, i);
        vector_of_global_cols.push_back(v);
    }*/

    std::vector<std::vector<int>> vector_of_global_cols(this->obs_mat_sys_rows);

    for (const auto& triplet : this->concatenated_triplets)
    {
        if (triplet.row() < this->obs_mat_sys_rows)
        {
            vector_of_global_cols[triplet.row()].push_back(triplet.col());
        }
    }

    std::cout << "vector_of_global_cols.size = :\n" << vector_of_global_cols.size() << std::endl;

    //std::cout << "here-1" << std::endl;

    /*for (auto& vec : vector_of_global_cols)
    {
        for (auto& num : vec)
        {
            std::cout << num << "\t";
        }
        std::cout << "\n";
    }*/

    //average

    //std::cout << "without files:\n";

    //for (size_t i = 0; i < this->amplitudes.size(); i++) // TODO : warning C4018: <: несоответствие типов со знаком и без знака
    for (Eigen::Index i = 0; i < this->amplitudes.size(); i++)
    {
        auto m = std::exp(this->amplitudes[i]);
        for (size_t j = 0; j < n_factors; j++)
        {
            /*std::cout << vector_of_global_cols[i][j] << " ";*/
            m /= this->solution[vector_of_global_cols[i][j]];
        }
        //std::cout << "ampl[i]: " << amplitudes[i] << "\tcorrection[i]: " << m << "\n";
        corrections.push_back((m));
    }

    //std::cout << "here-2\n";

    /*std::vector<size_t> vector_of_sizes =
    {
        this->sou_uniques->size(),
        this->rec_uniques->size(),
        this->cdp_uniques->size(),
        this->off_uniques->size(),
        this->az_uniques->size()
    };*/

    /*std::cout << "with files:\n";

    for (size_t i = 0; i < amplitudes.size(); i++)
    {
        double m = amplitudes[i];
        int shift{};
        for (size_t j = 0; j < n_factors; j++)
        {
            auto file_name = "solution_lsqr_" + files_names[j] + ".txt";
            m -= findByKey(file_name, vector_of_global_cols[i][j] - shift);
            shift += vector_of_sizes[j];
        }
        //std::cout << "ampl[i]: " << amplitudes[i] << "\tcorrection[i]: " << m << "\n";
    }*/

    //delete sou_uniques;
    this->sou_uniques.reset();
    //delete rec_uniques;
    this->rec_uniques.reset();
    //delete cdp_uniques;
    this->cdp_uniques.reset();
    //delete off_uniques;
    this->off_uniques.reset();
    //delete az_uniques;
    this->az_uniques.reset();

    auto end_entering_corrections = std::chrono::high_resolution_clock::now();

    auto duration_entering_corrections_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_entering_corrections - start_entering_corrections);

    std::cout << "\nCalculate corrections time: \t" << duration_entering_corrections_ms.count()
        << " milliseconds" << std::endl;

    /*=============================================================================================
                                 II. Writing corrections to output file
    =============================================================================================*/

    auto start_writing_corrections_to_output_file = std::chrono::high_resolution_clock::now();

    std::cout << "\nLOG:\tII. WRITING CORRECTIONS TO OUTPUT FILE\n";

    writeVectorToFile(corrections, "corrs.bin");

    auto end_writing_corrections_to_output_file = std::chrono::high_resolution_clock::now();

    auto duration_writing_corrections_to_output_file_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_writing_corrections_to_output_file - start_writing_corrections_to_output_file);

    std::cout << "\nWriting corrections to output file time: \t" << duration_writing_corrections_to_output_file_ms.count()
        << " milliseconds" << std::endl;
}

#endif // AMPLITUDES_CORRECTION_H