#include <iostream>
#include <queue>

#include "Prover.h"
#include "Verifier.h"

int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);
    // size_t matrix_size = 32;
    if (argc < 1) {
        std::cerr << "Usage: " << argv[0] << " <matrix_size> <filename>" << std::endl;
        return 1;
    }
    std::string filename = "";
    if (argc > 2) {
        // matrix_size = atoi(argv[1]);
        filename = argv[2];
    }

    std::queue<size_t> messageQueue;

    Prover prover(messageQueue, filename);
    Verifier verifier(messageQueue);

    uint8_t prover_status = prover.run();
    uint8_t verifier_status = verifier.run();

    while (prover_status != 2 &&
           verifier_status != 2 &&
           prover_status != 0 &&
           verifier_status != 0) {
        prover_status = prover.run();
        verifier_status = verifier.run();
    }

    return 0;
}