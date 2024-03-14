#include "Prover.h"

// #include <omp.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

void Prover::generateRandomBoolMatrix(size_t size) {
    this->matrix.resize(size, vector<bool>(size, false));
    random_device rd;
    mt19937 gen(rd());
    bernoulli_distribution dist(0.5);

    // #pragma omp parallel for
    for (size_t i = 0; i < size; ++i) {
        // #pragma omp parallel for
        for (size_t j = i; j < size; ++j) {
            if (i != j) {
                bool value = dist(gen);
                this->matrix[i][j] = value;
                this->matrix[j][i] = value;
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

        this->matrix[a][b] = this->matrix[b][a] = 1;
        this->matrix[a][c] = this->matrix[c][a] = 1;
        this->matrix[b][c] = this->matrix[c][b] = 1;
    }
}

vector<vector<bool>> Prover::generateCombinations(size_t n) {
    size_t totalCombinations = pow(2, n);
    vector<vector<bool>> result;

    for (size_t i = 0; i < totalCombinations; ++i) {
        vector<bool> combination;
        for (size_t j = 0; j < n; ++j) {
            // Shift i by j bits and check if the last bit is 1
            combination.push_back((i >> j) & 1);
        }
        reverse(combination.begin(), combination.end());
        result.push_back(combination);
    }

    return result;
}

void Prover::generateCombinations() {
    size_t n = log2(this->matrix.size()) * 2;
    size_t totalCombinations = pow(2, n);

    for (size_t i = 0; i < totalCombinations; ++i) {
        vector<bool> combination;
        for (size_t j = 0; j < n; ++j) {
            // Shift i by j bits and check if the last bit is 1
            combination.push_back((i >> j) & 1);
        }
        reverse(combination.begin(), combination.end());
        wss.push_back(combination);
    }
}

void Prover::generateCombinations(size_t n, const vector<size_t>& prevals,
                                  vector<vector<size_t>>& result) {
    //
    size_t totalCombinations = 1 << (n - prevals.size());
    result.resize(totalCombinations);
    for (size_t i = 0; i < totalCombinations; ++i) {
        result[i] = prevals;
        result[i].resize(n);

        for (size_t j = 0; j < n - prevals.size(); ++j) {
            result[i][prevals.size() + j] = (i & (1 << j)) != 0;
        }
        // result.push_back(combination);
    }

    // for (size_t i = 0; i < result.size(); ++i) {
    //     // vector<size_t> combination = prevals;
    //     // combination.resize(n);
    //     result[i].resize(n);
    //     for (size_t j = 0; j < n - prevals.size(); ++j) {
    //         result[i][prevals.size() + j] = (i & (1 << j)) != 0;
    //     }
    //     // result[i] = combination;
    // }
}

void Prover::generateGroupedCombinations(size_t n, const vector<size_t>& rs,
                                         size_t k, vector<vector<vector<size_t>>>& allCombinations) {
    // Expand A to size n

    vector<size_t> A(n, Zp);
    copy(rs.begin(), rs.end(), A.begin());
    uint8_t count = 0;
    uint8_t rs_pushed_count = 0;

    for (size_t i = 0; i < n; i += k) {
        vector<size_t> prevals{};
        vector<size_t> group(A.begin() + i, A.begin() + (i + k));
        for (const auto& val : group) {
            if (val != Zp) {
                prevals.push_back(val);
            }
        }

        // if (rs_pushed_count < rs.size()) {
        //     uint8_t end = min((i + k), n);
        //     auto ptr1 = rs.begin() + i;
        //     auto ptr2 = rs.begin() + end;
        //     rs_pushed_count += end;
        //     // copy values from rs to prevals
        //     copy(ptr1, ptr2, prevals.begin());
        // }

        // resize allCombinations[count]
        // size_t totalCombinations = 1 << (n - prevals.size());
        // allCombinations[count].resize(totalCombinations);

        generateCombinations(k, prevals, allCombinations[count++]);
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

    return this->matrix[a][b];
}

size_t Prover::poly_term(const vector<bool>& ws, const vector<size_t>& xs) {
    size_t res = 0;

    if (ws[0] == 0) {
        res = subZp(1, xs[0]);
    } else {
        res = xs[0];
    }

    for (size_t i = 1; i < ws.size(); ++i) {
        if (ws[i] == 0) {
            res = mulZp(res, subZp(1, xs[i]));
        } else {
            res = mulZp(res, xs[i]);
        }
    }

    return modZp(res);
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
    // vector<size_t>& rs = this->random_values;
    size_t k = log2(this->matrix.size());
    size_t n = k * 3;
    // #pragma omp parallel for
    for (size_t des_input = 0; des_input <= 2; des_input++) {
        // copy rs into temp vector and add des_input
        vector<size_t> temp_rs(random_values.begin(), random_values.end());
        temp_rs.push_back(des_input);

        vector<vector<vector<size_t>>> groupedCombinations{vector<vector<size_t>>{},
                                                           vector<vector<size_t>>{},
                                                           vector<vector<size_t>>{}};

        generateGroupedCombinations(n, temp_rs, k, groupedCombinations);

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

    size_t k = log2(this->matrix.size());
    size_t n = k * 3;

    vector<vector<vector<size_t>>> groupedCombinations{vector<vector<size_t>>{},
                                                       vector<vector<size_t>>{},
                                                       vector<vector<size_t>>{}};
    generateGroupedCombinations(n, random_values, k, groupedCombinations);

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
        return matrix;  // return empty matrix on failure
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
    cout << " ========= prover turn ============ " << endl;
    if (this->round == 0) {
        // first round
        // calculate the number of triangles
        size_t triangles = count_triangles();
        cout << "Prover: Number of triangles: " << triangles << endl;

        // send the number of triangles to the verifier
        this->messageQueue.push(triangles);
        this->round++;
        return 1;
    } else if (this->round == 1) {
        // second round
        // calculate the three values vector
        vector<size_t> subPolyMle = calculate_subPolyMle();
        cout << "Prover in "
             << "round " << round << " subPolyMle: "
             << subPolyMle[0] << " "
             << subPolyMle[1] << " "
             << subPolyMle[2] << endl;

        // send the three values to the verifier
        for (const size_t& val : subPolyMle) {
            this->messageQueue.push(val);
        }
        this->round++;
        return 1;
    }

    // other rounds
    // get the random value from verifier
    // cout << "where ? 1" << endl;
    size_t r = this->messageQueue.front();
    // cout << "where ? 2" << endl;
    this->messageQueue.pop();
    // cout << "where ? 3" << endl;

    cout << "Prover in "
         << " round " << this->round
         << " received r: " << r
         << endl;

    // push to random_values
    this->random_values.push_back(r);

    // calculate the subPolyMle
    cout << "here again?" << endl;
    vector<size_t> subPolyMle = calculate_subPolyMle();
    cout << "here again?2" << endl;

    cout << "Prover in "
         << "round " << this->round << " subPolyMle: "
         << subPolyMle[0] << " "
         << subPolyMle[1] << " "
         << subPolyMle[2] << endl;

    // send the three values to the verifier
    for (const size_t& val : subPolyMle) {
        this->messageQueue.push(val);
    }

    this->round++;
    if (this->random_values.size() == this->xs_size) {
        cout << "Prover: Done" << endl;
        return 2;
    }

    return 1;
}

// if (i < rs.size()) {
//     for (const auto& r : rs) {
//         prevals.push_back(r);
//         if (prevals.size() == k) {
//             break;
//         }
//     }
// }

// vector<size_t> prevals{};
// if (rs.size() == 0) {
// } else if (i + k < rs.size()) {
//     copy(rs.begin() + i, rs.begin() + (i + k), prevals.begin());
// } else if (i < rs.size()) {
//     // prevals = vector<size_t>(rs.begin() + i, rs.end() - 1);
//     copy(rs.begin() + i, rs.end() - 1, prevals.begin());
// }