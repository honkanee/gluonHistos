#include "GluonHistosFill.h"

#include <ctime>
#include <chrono>

int main() {

    auto start = chrono::system_clock::now(); //timer

    GluonHistosFill g;

    g.Loop();

    // timer
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_sec = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "finished at " << std::ctime(&end_time)
            << "elapsed time: " << elapsed_sec.count()/60 << "min" << endl;
};
