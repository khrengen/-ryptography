#include <iostream>
#include <vector>
#include <stack>
#include <ctime>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmpxx.h>

#include "qs.cpp"

const static uint32_t TRIAL_THRESHOLD = 1000000000;

// Maximum input size we can handle (bits).
const static uint32_t MAX_DIGITS = 100;


int main() {

    mpz_class N;
    while (std::cin >> N) {
        if (mpz_sizeinbase(N.get_mpz_t(), 2) > MAX_DIGITS) {
            std::cout << "fail" << std::endl << std::endl; // Too many digits.
            continue;
        }

        if (mpz_probab_prime_p(N.get_mpz_t(), 10)) {
            // N is prime.
            continue;
        }

        std::stack<mpz_class> factors;
        factors.push(N);

        while (!factors.empty()) {
            mpz_class factor = factors.top();
            factors.pop();

            if (mpz_probab_prime_p(factor.get_mpz_t(), 10)) {
                // N is prime.
                std::cout << factor << std::endl;
                continue;
            }

            
                mpz_class result = quadraticSieve(factor);
                factors.push(result);
                factors.push(factor/result);
            
        }
        std::cout << std::endl;
    }
    return 0;
}

//-lgmpxx -lgmp to compile
