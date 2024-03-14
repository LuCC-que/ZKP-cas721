#include "Verifier.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

uint8_t Verifier::run() {
    cout << " ========= verifier turn ============ " << endl;
    // check if there are three values or one value in message Queue
    if (this->messageQueue.size() == 1) {
        // push this value into the stack
        size_t c1 = this->messageQueue.front();
        this->verify_stack.push(c1);
        cout << "verifier received a value: " << c1 << endl;
        this->messageQueue.pop();

        return 1;
    } else if (messageQueue.size() < 3) {
        cout << "Error: messageQueue.size() < 3" << endl;
        return 0;
    }

    size_t c1 = this->messageQueue.front();
    this->messageQueue.pop();
    size_t c2 = this->messageQueue.front();
    this->messageQueue.pop();
    size_t c3 = this->messageQueue.front();
    this->messageQueue.pop();

    // check if verify_stack is not empty
    if (verify_stack.empty()) {
        cout << "Error: verify_stack is empty" << endl;
        return 0;
    }

    // if c1 + c2 is the value in verify_stack
    size_t verify_stack_top = verify_stack.top();
    verify_stack.pop();

    if (c1 == c2 && c2 == c3) {
        if (c1 == verify_stack_top) {
            cout << "verified c1 == c2 == c3 == verify_stack.top()"
                 << endl
                 << "c1, c2, c3: " << c1
                 << " verify_stack.top(): "
                 << verify_stack_top
                 << endl;

            cout << "Verifier: Done" << endl;

            return 2;
        } else {
            cout << "Error: c1 == c2 == c3 != verify_stack.top()"
                 << endl
                 << "c1, c2, c3: " << c1
                 << " verify_stack.top(): "
                 << verify_stack_top
                 << endl;
            return 0;
        }
    }

    size_t c1_c2 = addZp(c1, c2);
    if (c1_c2 != verify_stack_top) {
        // push c3 into verify_stack
        cout << "Error: c1 + c2 != verify_stack.top()"
             << endl
             << "c1 + c2: " << c1_c2
             << " verify_stack.top(): "
             << verify_stack_top
             << endl;
        return 0;
    }

    cout << "verified c1 + c2 == verify_stack.top()"
         << endl
         << "c1 + c2: " << c1_c2
         << " verify_stack.top(): "
         << verify_stack_top
         << endl;

    // pick a random value from 1 to Zp
    static random_device rd;
    static mt19937 gen(rd());
    static uniform_int_distribution<int> distribution(1, Zp);
    size_t r = distribution(gen);

    // calculate the interpolation
    size_t result = interpolation(c1, c2, c3, r);

    cout << "interpolation result: " << result << endl;

    // save the result into verify_stack
    this->verify_stack.push(result);

    // push the r1 into messageQueue
    this->messageQueue.push(r);

    return 1;
}