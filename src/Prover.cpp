#include "Prover.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

void Prover::generateRandomBoolMatrix(size_t size) {
    matrix.resize(size, vector<bool>(size, false));
    random_device rd;
    mt19937 gen(rd());
    bernoulli_distribution dist(0.5);

    // #pragma omp parallel for
    for (size_t i = 0; i < size; ++i) {
        // #pragma omp parallel for
        for (size_t j = i; j < size; ++j) {
            if (i != j) {
                bool value = dist(gen);
                matrix[i][j] = value;
                matrix[j][i] = value;
            }
        }
    }

    // Ensure at least one triangle
    if (size >= 3) {
        vector<int> indices(size);
        iota(indices.begin(), indices.end(), 0);
        shuffle(indices.begin(), indices.end(), gen);

        // Select the first three indices to form a triangle
        int a = indices[0];
        int b = indices[1];
        int c = indices[2];

        matrix[a][b] = matrix[b][a] = 1;
        matrix[a][c] = matrix[c][a] = 1;
        matrix[b][c] = matrix[c][b] = 1;
    }

    cout << "Generated matrix:" << endl;
    for (const auto& row : matrix) {
        for (bool elem : row) {
            cout << (elem ? "1" : "0") << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void Prover::generateCombinations() {
    /*
        Explaining (i & (1 << j)) != 0
        the wss.size() is total the number of combinations
        the wss[0].size() is the number of variables

        for example if wss[0].size() = 3
        then the combinations are 2^3 = 8

        000
        001
        010
        ...
        111

        our goal is to generate all the combinations
        but in vector<bool> form, we can acheive this
        by shift least significant bit of value 1 by j
        position and then do bitwise AND with i

        like this, suppose size = 3

        i = 1, we want vector<bool> = {0, 0, 1}


        i: 1, j: 0, 1 << j: 1, i & (1 << j): 1, result: 1
        i: 1, j: 1, 1 << j: 2, i & (1 << j): 0, result: 0
        i: 1, j: 2, 1 << j: 4, i & (1 << j): 0, result: 0

        more intuitive example

        1 = 001 in 3 bits
        001 & 001 = 1
        001 & 010 = 0
        001 & 100 = 0

        rest of the case are like this

    */
    for (size_t i = 0; i < wss.size(); i++) {
        for (size_t j = 0; j < wss[0].size(); j += 4) {
            // I don't have support for 2x2 matrix
            // so I am assuming the size is at least 4
            // and multiple of 4, so I can do the
            // loop unrolling to acclelerate the process
            wss[i][j] = (i & (1 << j)) != 0;
            wss[i][j + 1] = (i & (1 << (j + 1))) != 0;
            wss[i][j + 2] = (i & (1 << (j + 2))) != 0;
            wss[i][j + 3] = (i & (1 << (j + 3))) != 0;
        }
    }
}

void Prover::generateCombinations(size_t n, const vector<size_t>& prevals,
                                  uint8_t count) {
    /*
        this function is similar to the above function

        but here we have fixed values

        like if we have size = 3, but prevals = {124}

        then we have to generate the combinations of 0s and 1s
        of size 3 - 1 = 2

        {124, 0, 0}
        {124, 0, 1}
        {124, 1, 0}
        {124, 1, 1}

        and this will be the all permutations of a group
    */
    size_t totalCombinations = 1 << (n - prevals.size());
    groupedCombinations[count].resize(totalCombinations);
    for (size_t i = 0; i < totalCombinations; ++i) {
        groupedCombinations[count][i].resize(n);

        // copy the fixed values into groupedCombination
        // to fixed the unchanged position
        copy(prevals.begin(), prevals.end(), groupedCombinations[count][i].begin());
        for (size_t j = 0; j < n - prevals.size(); ++j) {
            groupedCombinations[count][i][prevals.size() + j] = (i & (1 << j)) != 0;
        }
    }
}

void Prover::generateGroupedCombinations(size_t n, const vector<size_t>& rs,
                                         size_t k) {
    /*
    We separte the random values into k groups

    since we are always have three group

    take 4x4 as example
    k = log2(4)*2, n = log2(4)*3, put it simply
    log2(4) is number of elements in a group
    k is the number of elements in a group combination
    n is the number of elements of three vertices

    x1x2x3x4
    x1x2x5x6
    x3x4x5x6

    xi is variable that can be 0 or 1,

    thus, we can have x1x2 as G1, x3x4 as G2, x5x6 as G3

    xi can be replaced by fixed value from received random values, then it will
    not have chance to be 0 or 1, it just fixed value

    k = log2(matrix.size())*2, represent how much many elements
    in a group, in 8x8, k = 6, n = 9

    x1x2x3x4x5x6
    x1x2x3x7x8x9
    x4x5x6x7x8x9

    this works for any size of matrix, 64x64, 128x128, etc
    but just slower for larger matrix

    */

    uint8_t count = 0;
    uint8_t rs_pushed_count = 0;
    const size_t rs_size = rs.size();

    for (size_t i = 0; i < n; i += k) {
        vector<size_t> prevals{};

        // push the random values into prevals as fixed value
        //  if have position left that will be free values that can be 0 or 1
        for (; rs_pushed_count < rs_size;) {
            prevals.push_back(rs[rs_pushed_count++]);
            size_t prevals_size = prevals.size();
            if (prevals_size == k) {
                break;
            }
        }

        generateCombinations(k, prevals, count++);
    }
}

bool Prover::f(const vector<bool>& xs) {
    size_t a = 0;
    size_t b = 0;
    size_t fNv = xs.size();

    for (size_t i = fNv, j = fNv / 2, div = 1;
         i > (fNv / 2) && j > 0; --i, --j, div *= 2) {
        a += div * xs[j - 1];
        b += div * xs[i - 1];
    }

    return matrix[a][b];
}

size_t Prover::poly_term(const vector<bool>& ws, const vector<size_t>& xs) {
    size_t res = 1;

    for (size_t i = 0; i < ws.size(); i++) {
        res = mulZp(res, !ws[i] * subZp(1, xs[i]) + ws[i] * xs[i]);
    }

    return res;
}

size_t Prover::F(const vector<size_t>& xs) {
    size_t res = 0;

    for (const vector<bool>& ws : wss) {
        res = addZp(res, mulZp(f(ws), poly_term(ws, xs)));
    }

    return modZp(res);
}

vector<size_t> Prover::calculate_subPolyMle() {
    vector<size_t> count = {0, 0, 0};
    // vector<size_t>& rs = random_values;
    size_t k = log2(matrix.size());
    size_t n = k * 3;

    for (size_t des_input = 0; des_input <= 2; des_input++) {
        // copy rs into temp vector and add des_input
        vector<size_t> temp_rs(random_values.begin(), random_values.end());
        temp_rs.push_back(des_input);

        /*
         Chek the description of generateGroupedCombinations

         since we have three groups, we have to generate the combinations

         take 4x4 as example

         x1x2 = G1 = {all permutation of G1}
         x3x4 = G2 = {all permutation of G2}
         x5x6 = G3 = {all permutation of G3}

         Then we can pick one possibilty from each G1, G2, G3
         denote g1, g2, g3

         then we can concanate g1g2, g1g3, g2g3 to be all the possible
         inputs for this round, and so on


        */

        generateGroupedCombinations(n, temp_rs, k);

        for (const auto& group1 : groupedCombinations[0]) {
            for (const auto& group2 : groupedCombinations[1]) {
                for (const auto& group3 : groupedCombinations[2]) {
                    vector<size_t> v1(group1.begin(), group1.end());
                    v1.insert(v1.end(), group2.begin(), group2.end());
                    vector<size_t> v2(group1.begin(), group1.end());
                    v2.insert(v2.end(), group3.begin(), group3.end());
                    vector<size_t> v3(group2.begin(), group2.end());
                    v3.insert(v3.end(), group3.begin(), group3.end());

                    count[des_input] = addZp(count[des_input],
                                             mulZp(mulZp(F(v1), F(v2)), F(v3)));
                }
            }
        }
    }

    return count;
}

size_t Prover::count_triangles() {
    size_t count = 0;

    size_t k = log2(matrix.size());
    size_t n = k * 3;

    generateGroupedCombinations(n, random_values, k);

    for (const auto& group1 : groupedCombinations[0]) {
        for (const auto& group2 : groupedCombinations[1]) {
            for (const auto& group3 : groupedCombinations[2]) {
                vector<size_t> v1(group1.begin(), group1.end());
                v1.insert(v1.end(), group2.begin(), group2.end());
                vector<size_t> v2(group1.begin(), group1.end());
                v2.insert(v2.end(), group3.begin(), group3.end());
                vector<size_t> v3(group2.begin(), group2.end());
                v3.insert(v3.end(), group3.begin(), group3.end());

                count = addZp(count, mulZp(mulZp(F(v1), F(v2)), F(v3)));
            }
        }
    }

    return count;
}

vector<vector<bool>> Prover::readMatrixFromFile(const string& filename) {
    vector<vector<bool>> matrix;
    ifstream inFile(filename);
    string line;

    if (!inFile.is_open()) {
        cerr << "Failed to open " << filename << endl;
        return matrix;
    }

    while (getline(inFile, line)) {
        istringstream iss(line);
        vector<bool> row;
        bool val;

        while (iss >> val) {
            row.push_back(val);
        }

        matrix.push_back(row);
    }

    inFile.close();
    return matrix;
}

uint8_t Prover::run() {
    cout << "In round " << round << endl;
    if (round == 0) {
        // first round
        // calculate the number of triangles
        size_t triangles = count_triangles();
        cout << "P: Initial G " << triangles << endl;

        // send the number of triangles to the verifier
        cout << "P: sends " << triangles << " to V" << endl;
        messageQueue.push(triangles);
        round++;
        return 1;
    } else if (round == 1) {
        // second round
        // calculate the three values vector
        vector<size_t> subPolyMle = calculate_subPolyMle();

        cout << "P: sends "
             << subPolyMle[0] << " "
             << subPolyMle[1] << " "
             << subPolyMle[2] << endl;

        // send the three values to the verifier
        for (const size_t& val : subPolyMle) {
            messageQueue.push(val);
        }
        round++;
        return 1;
    }

    // other rounds
    // get the random value from verifier
    size_t r = messageQueue.front();
    messageQueue.pop();

    cout << "P:"
         << " received r " << r
         << " in round "
         << endl;

    random_values.push_back(r);

    vector<size_t> subPolyMle = calculate_subPolyMle();

    cout << "P: sends "
         << subPolyMle[0] << " "
         << subPolyMle[1] << " "
         << subPolyMle[2] << " to V in round "
         << endl;

    // send the three values to the verifier
    for (const size_t& val : subPolyMle) {
        messageQueue.push(val);
    }

    round++;
    if (random_values.size() == xs_size) {
        cout << "P: Done sending" << endl;
        return 2;
    }

    return 1;
}
