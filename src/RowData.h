//
// Created by Gregory Johnson on 6/21/24.
//

#ifndef ROWDATA_H
#define ROWDATA_H
#include<string>
#include<vector>

struct RowValues {
    RowValues(const std::string &name, const float value){
        this->name = name;
        this->value = value;
    }
    std::string name;
    float value;
};
struct RowData {
    std::string name;
    std::vector<float> rowValues;
    static bool SortComparison(const RowData& a, const RowData& b) {
        return std::stoi(a.name) < std::stoi(b.name);
    }
};

#endif //ROWDATA_H
