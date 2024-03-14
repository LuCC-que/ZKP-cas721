#include <chrono>
#include <iostream>
#include <memory>
#include <queue>

#include "Prover.h"
#include "Verifier.h"

using namespace std::chrono;

int main(int argc, char* argv[]) {
    auto start = high_resolution_clock::now();
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    size_t matrix_size = 4;
    if (argc < 2) {
        std::cerr << std::endl
                  << "Usage: " << argv[0] << " <matrix_size> <filename>"
                  << std::endl;

        std::cerr << "if you have filename, matrix_size will be ignored"
                  << std::endl;

        std::cerr << "I also did not build support for 2x2 matrix"
                  << std::endl
                  << "so please use matrix_size >= 4"
                  << std::endl
                  << std::endl;
        return 1;
    }
    std::string filename = "";

    if (argc == 2) {
        matrix_size = atoi(argv[1]);
    }
    if (argc > 2) {
        // matrix_size = atoi(argv[1]);
        filename = argv[2];
    }

    std::queue<size_t> messageQueue;

    std::unique_ptr<Prover> prover;
    Verifier verifier(messageQueue);

    // Dynamically allocate the Prover based on the condition
    if (filename == "") {
        prover = std::make_unique<Prover>(messageQueue, matrix_size);
    } else {
        prover = std::make_unique<Prover>(messageQueue, filename);
    }

    // Now, use the prover and verifier
    uint8_t prover_status = prover->run();
    uint8_t verifier_status = verifier.run();

    while (prover_status != 2 && verifier_status != 2 &&
           prover_status != 0 && verifier_status != 0) {
        prover_status = prover->run();
        verifier_status = verifier.run();
    }

    auto stop = high_resolution_clock::now();

    // Calculate duration
    auto duration = duration_cast<microseconds>(stop - start);

    cout << endl
         << "Time taken by protocol: "
         << duration.count()
         << " microseconds"
         << endl;

    return 0;
}