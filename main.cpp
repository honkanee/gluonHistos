#include "GluonHistosFill.h"
#include "drawGluonHistos.C"

int main() {
    GluonHistosFill g;

    g.Loop();
    plot();
};
