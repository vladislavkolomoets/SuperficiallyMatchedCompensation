#ifndef CORRECTPROC_H
#define CORRECTPROC_H

#include <unordered_map>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

std::unordered_map<int, double> readFromFileAsDict(const std::string& _file_name)
{
    std::unordered_map<int, double> values{};
    std::ifstream file(_file_name);
    std::string line{};

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        int a{};
        double b{};
        char comma{};

        if (!(iss >> a >> comma >> b))
        {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue;
        }
        values[a] = b;
    }

    return values;
}

int findByKey(const std::string& _file_name, int key)
{
    std::unordered_map<int, double> valueMap = readFromFileAsDict(_file_name);
    auto it = valueMap.find(key);
    if (it != valueMap.end())
    {
        return it->second; // TODO : warning C4244: return: преобразование "_Ty2" в "int", возможна потеря данных
    }
    else
    {
        //std::cerr << "Value for a key (" << key << ") not found in the file: " << _file_name << "\n";
        return 0;
    }
}

#endif // LIBRARY_H