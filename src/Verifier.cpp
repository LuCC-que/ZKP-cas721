#include "Verifier.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

uint8_t Verifier::run() {
    // check if there are three values or one value in message Queue
    if (messageQueue.size() == 1) {
        // push this value into the stack
        size_t c1 = messageQueue.front();
        verify_stack.push(c1);
        cout << "V: next verified value "
             << c1
             << endl
             << endl;
        messageQueue.pop();

        return 1;
    } else if (messageQueue.size() < 3) {
        cout << "Error: messageQueue.size() < 3" << endl;
        return 0;
    }

    size_t c1 = messageQueue.front();
    messageQueue.pop();
    size_t c2 = messageQueue.front();
    messageQueue.pop();
    size_t c3 = messageQueue.front();
    messageQueue.pop();

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
            cout << "V: "
                 << "c1, c2, c3 " << c1
                 << ", last verified value "
                 << verify_stack_top
                 << " ------ verify accepted"
                 << endl;

            cout << "Verifier: Proof accepted" << endl;

            return 2;
        } else {
            cout << "Error: c1 == c2 == c3 != verify_stack.top()"
                 << endl
                 << "c1, c2, c3 " << c1
                 << " last verified value"
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
             << "c1 + c2 " << c1_c2
             << ", last verified value "
             << verify_stack_top
             << endl;
        return 0;
    }

    cout << "V: "
         << "(c1 + c2) mod Zp " << c1_c2
         << ", last verified value "
         << verify_stack_top
         << " ------ verify accepted"
         << endl;

    // pick a random value from 1 to Zp
    static random_device rd;
    static mt19937 gen(rd());
    static uniform_int_distribution<int> distribution(1, Zp);
    size_t r = distribution(gen);

    // calculate the interpolation
    size_t result = interpolation(c1, c2, c3, r);

    cout << "V: next verified value  " << result << endl
         << endl;

    // save the result into verify_stack
    verify_stack.push(result);

    // push the r1 into messageQueue
    messageQueue.push(r);

    return 1;
}