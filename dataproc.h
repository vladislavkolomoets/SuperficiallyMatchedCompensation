#ifndef DATAPROC_H
#define DATAPROC_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <unordered_map>
//#include <chrono>
#include <map>
#include <memory>
#include <algorithm>

using EnumerationNumber = std::map<int, size_t>;


std::unique_ptr<int[]> readColomnFromBinaryFile(const std::string& _file_name, const size_t _rows_size)
{
    std::ifstream file(_file_name, std::ios::in | std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "error with openning file: " << _file_name << std::endl;
        return {};
    }

    //int* data = new int[_rows_size];
    auto data = std::make_unique<int[]>(_rows_size);
    //file.read((char*)data, sizeof(int) * _rows_size);
    file.read(reinterpret_cast<char*>(data.get()), sizeof(int) * _rows_size);

    file.close();

    return data;
}

EnumerationNumber calcEnumerationNumber(const int* _coordinates, const size_t _size)
{
    // Valid check
    if (_coordinates == NULL)
    {
        return{};
    }

    EnumerationNumber coordinates_map{};

    size_t current_index = 0;

    // If the given coordinate is already in the map, we use its ordinal number 
    // Otherwise, a new serial number is assigned 

    for (size_t i = 0; i < _size; i++)
    {
        //const CoordinateInfo& coord = coordinate;
        auto it = coordinates_map.find(_coordinates[i]);
        if (it == coordinates_map.end())
        {
            coordinates_map[_coordinates[i]] = current_index;
            current_index++;
        }

    }

    return coordinates_map;
}

std::vector<size_t> makeCodes(const EnumerationNumber& _coords_map, const int* _coords, const size_t _rows_size)
{
    //std::cout << "LOG\nIn makeCodes()\n";

    if (_coords_map.size() == 0 || _coords == NULL)
    {
        return {};
    }
    auto number_of_rows = _rows_size;

    std::vector<size_t> codes{};
    for (size_t i = 0; i < number_of_rows; i++)
    {
        auto it = _coords_map.find(_coords[i]);
        codes.emplace_back(it->second);
    }

    return codes;
}


void writeVectorToFile(const std::vector<double>& _vector, const std::string& _file_name)
{
    if (_vector.empty())
    {
        return;
    }

    std::ofstream file(_file_name, std::ios::out | std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "error with opening file: " << _file_name << std::endl;
        return;
    }

    file.write(reinterpret_cast<const char*>(_vector.data()), _vector.size() * sizeof(double));

    file.close();

    std::cout << "\nVector successfully written to file: " << _file_name << std::endl;
}



#endif // LIBRARY_H