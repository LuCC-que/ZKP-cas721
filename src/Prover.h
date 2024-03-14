#ifndef PROVER_H
#define PROVER_H

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <vector>

#define Zp 524287
#define modZp(x) (((x) % Zp))
#define addZp(x, y) (modZp(modZp(x) + modZp(y)))
#define subZp(x, y) (modZp((modZp(x) + Zp) - modZp(y)))
#define mulZp(x, y) (modZp(modZp(x) * modZp(y)))

using namespace std;
class Prover {
   private:
    vector<vector<bool>> matrix;
    vector<size_t> random_values;
    queue<size_t>& messageQueue;

    vector<vector<bool>> wss;

    size_t matrix_size;
    size_t xs_size;
    size_t round = 0;

    void generateRandomBoolMatrix(size_t size);
    inline vector<vector<bool>> generateCombinations(size_t n);
    void generateCombinations();
    inline void generateCombinations(size_t n, const vector<size_t>& prevals,
                                     vector<vector<size_t>>& result);
    void generateGroupedCombinations(size_t n, const vector<size_t>& rs, size_t k,
                                     vector<vector<vector<size_t>>>& allCombinations);

    inline bool f(const vector<bool>& xs);
    inline size_t poly_term(const vector<bool>& ws, const vector<size_t>& xs);

    inline size_t F(const vector<size_t>& xs);

    inline vector<size_t> calculate_subPolyMle();
    inline size_t count_triangles();

    vector<vector<bool>> readMatrixFromFile(const string& filename);

   public:
    Prover(queue<size_t>& messageQueue, size_t matrix_size) : messageQueue(messageQueue), matrix_size(matrix_size) {
        generateRandomBoolMatrix(matrix_size);
        xs_size = log2(matrix_size) * 3;
        generateCombinations();
    }
    Prover(queue<size_t>& messageQueue, string matrixFile) : messageQueue(messageQueue) {
        matrix = readMatrixFromFile(matrixFile);
        matrix_size = matrix.size();
        xs_size = log2(matrix_size) * 3;
        generateCombinations();
    }
    ~Prover() {}

    uint8_t run();
};

#endif