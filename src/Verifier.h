#ifndef VERIFIER_H
#define VERIFIER_H

#include <cstdint>
#include <queue>
#include <stack>
#include <vector>

#define Zp 524287
#define modZp(x) (((x) % Zp))
#define addZp(x, y) (modZp(modZp(x) + modZp(y)))
#define subZp(x, y) (modZp((modZp(x) + Zp) - modZp(y)))
#define mulZp(x, y) (modZp(modZp(x) * modZp(y)))

#define first_term(x, y) (mulZp(mulZp(mulZp(262144, x), subZp(y, 1)), subZp(y, 2)))
#define second_term(x, y) (mulZp(mulZp(x, y), subZp(y, 2)))
#define thrid_term(x, y) (mulZp(mulZp(mulZp(262144, x), y), subZp(y, 1)))

#define interpolation(c1, c2, c3, r) (addZp(subZp(first_term(c1, r),   \
                                                  second_term(c2, r)), \
                                            thrid_term(c3, r)))

using namespace std;

class Verifier {
   private:
    queue<size_t>& messageQueue;
    vector<size_t> random_values;
    stack<size_t> verify_stack;

   public:
    Verifier(queue<size_t>& messageQueue) : messageQueue(messageQueue) {}
    ~Verifier() {}

    uint8_t run();
};

#endif