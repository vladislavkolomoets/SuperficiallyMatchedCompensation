#include "AmplitudesCorrection.h"
//#include "Deconvolution.h"

int main()
{
    //const size_t SIZE_SCAC_fin = 6147792; 

    std::pair<size_t, size_t> sizes_of_ampl_matrix = { ROWS_SIZE, 500 };

    std::vector<std::string> file_names =
    {
        "/mnt/c/Users/geouser7/Documents/SCAC_fin/data_SCAC_fin.bin",
        "/mnt/c/Users/geouser7/Documents/SCAC_fin/SN_SCAC_fin.bin",
        "/mnt/c/Users/geouser7/Documents/SCAC_fin/RN_SCAC_fin.bin",
        "/mnt/c/Users/geouser7/Documents/SCAC_fin/CDPN_SCAC_fin.bin",
        "/mnt/c/Users/geouser7/Documents/SCAC_fin/OFFN_SCAC_fin.bin",
        "/mnt/c/Users/geouser7/Documents/SCAC_fin/AZN_SCAC_fin.bin"
    };

    std::vector<std::string> file_names_for_VS =
    {
        "C:\\Users\\geouser7\\Documents\\SCAC_fin\\data_SCAC_fin.bin",
        "C:\\Users\\geouser7\\Documents\\SCAC_fin\\SN_SCAC_fin.bin",
        "C:\\Users\\geouser7\\Documents\\SCAC_fin\\RN_SCAC_fin.bin",
        "C:\\Users\\geouser7\\Documents\\SCAC_fin\\CDPN_SCAC_fin.bin",
        "C:\\Users\\geouser7\\Documents\\SCAC_fin\\OFFN_SCAC_fin.bin",
        "C:\\Users\\geouser7\\Documents\\SCAC_fin\\AZN_SCAC_fin.bin"
    };

    SuperficiallyMatchedCompensation<float> algorithm; // TODO : message : выполняется компиляция ссылки на экземпляр класс шаблон функции "SuperficiallyMatchedCompensation<float>"

    algorithm.processInputData(file_names, sizes_of_ampl_matrix, 4, 500, 1800, "RMSoutput.bin");

    algorithm.calculateCorrections();

    algorithm.enterCorrections(); // TODO : message : см. первую ссылку на "SuperficiallyMatchedCompensation<float>::enterCorrections" в "main"

    return 0;
}