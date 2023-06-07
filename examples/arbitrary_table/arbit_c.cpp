#include <iostream>
#include <fstream>

const int indepCount = 4; // independent variables count
const int depCount = 3;   // dependent variables count
const int d1 = 5, d2 = 9, d3 = 17, d4 = 33; // independent dimensions

int main() {
    int dims = d1 * d2 * d3 * d4;
    float table [depCount][dims];  // note the dimensions are reversed compared to fortran
    int k; // index counter
    int i1, i2, i3, i4;
    int j; // dependent variables

    for (j = 0; j < depCount; j++) {
        k = 0;
        for (i1 = 1; i1 <= d1; i1++) { // slowest indep
            for (i2 = 1; i2 <= d2; i2++) {
                for (i3 = 1; i3 <= d3; i3++) {
                    for (i4 = 1; i4 <= d4; i4++) { // fastest indep                   
                        table [j][k++] = (j+1) + 10 * i1 + 100 * i2 + 1000 * i3 + 10000 * i4; // a fake formula
                    }
                }
            }
        }
    }
    std::ofstream file("kookup_3dep_4indep_cpp.bin", std::ios::binary);
    if (!file) {
        std::cout << "Error: failed to open file to write" << std::endl;
        return 1;
    }
    std::cout<<sizeof(table)<<"\n";
    file.write(reinterpret_cast<char*>(table), sizeof(table));
    file.close();

    return 0;
}
