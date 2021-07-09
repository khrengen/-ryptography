#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>
 

#include "matrix.h"

// Minimal smoothness bound.
const static uint32_t MINIMAL_BOUND = 300;

// Sieving interval length.
const static uint32_t INTERVAL_LENGTH = 655360;


std::vector<uint32_t> generateFactorBase(const mpz_class& N, uint32_t B) {
    std::vector<uint32_t> factorBase;

    std::vector<bool> sieve(B + 1, false);
    for (uint32_t p = 2; p <= B; ++p) {
        if (sieve[p])
            continue;

        if (mpz_legendre(N.get_mpz_t(), mpz_class(p).get_mpz_t()) == 1)
            factorBase.push_back(p);

        for (uint32_t i = p; i <= B; i += p)
            sieve[i] = true;
    }

    return factorBase;
}

/*
 * Returns b^e (mod m) using right-to-left binary method.
 */
uint64_t modularPow(uint64_t b, uint64_t e, uint64_t m) {
    uint64_t result = 1;
    while (e > 0) {
        if (e & 1) 
            result = (result * b) % m; 
        e >>= 1;
        b = (b * b) % m; 
    }
    return result;
}

/*
 * Calculate the Legendre symbol.
 */
int32_t legendreSymbol(uint32_t a, uint32_t p) {
    uint64_t result = modularPow(a, (p - 1) / 2, p);
    return result > 1 ? -1 : result;
}


std::pair<uint32_t, uint32_t> tonelliShanks(uint32_t n, uint32_t p) {
    if (p == 2)
        return std::make_pair(n, n); // Double root.

    // Define Q2^S = p - 1.
    uint64_t Q = p - 1;
    uint64_t S = 0;
    while (Q % 2 == 0) {
        Q /= 2;
        ++S;
    }

    uint64_t z = 2;
    while (legendreSymbol(z, p) != -1)
        ++z;

    uint64_t c = modularPow(z, Q, p);            // c = z^Q         (mod p)
    uint64_t R = modularPow(n, (Q + 1) / 2, p);  // R = n^((Q+1)/2) (mod p)
    unsigned long long t = modularPow(n, Q, p);            // t = n^Q         (mod p)
    uint64_t M = S;

    while (t % p != 1) {
     
        int32_t i = 1;
        while (modularPow(t, std::pow(2, i), p) != 1) {
            ++i;
        }
        if (i >= M) {
            break;
        }

        uint64_t b = modularPow(c, std::pow(2, M - i - 1), p);
        R = (R * b) % p;    // R = Rb (mod p)
        t = (t * b * b) % p;  // t = tb^2
        c = (b * b) % p;      // c = b^2 (mod p)
        M = i;
    }

    return std::make_pair(R, p - R);
    

}

/*
 * A basic implementation of the Quadratic Sieve algorithm.
 */
mpz_class quadraticSieve(const mpz_class& N) {

    const float logN = mpz_sizeinbase(N.get_mpz_t(), 2) * std::log(2);
    const float loglogN = std::log(logN);
    const mpz_class sqrtN = sqrt(N);

    // Smoothness bound B.
    const uint32_t B = MINIMAL_BOUND + std::ceil(std::exp(0.55*std::sqrt(logN * loglogN)));

    /*
     * Step 1
     *
     * Generate factor base.
     */  
    const std::vector<uint32_t> factorBase = generateFactorBase(N, B);
    

    /*
     * Step 2
     *
     * Calculate start indices for each number in the factor base.
     */

    std::pair<std::vector<uint32_t>, std::vector<uint32_t> > startIndex(
        std::vector<uint32_t>(factorBase.size()), // Vector of first start index.
        std::vector<uint32_t>(factorBase.size())  // Vector of second start index.
    );
    for (uint32_t i = 0; i < factorBase.size(); ++i) {
        uint32_t p = factorBase[i];
        uint32_t N_mod_p = mpz_class(N % p).get_ui(); 

        std::pair<uint32_t, uint32_t> x = tonelliShanks(N_mod_p, p);
        startIndex.first[i] = mpz_class((((x.first - sqrtN) % p) + p) % p).get_ui();
        startIndex.second[i] = mpz_class((((x.second - sqrtN) % p) + p) % p).get_ui();
    }


    uint32_t intervalStart = 0;
    uint32_t intervalEnd = INTERVAL_LENGTH;

    std::vector<uint32_t> smooth;                     
    std::vector<std::vector<uint32_t> > smoothFactors; 
    std::vector<float> logApprox(INTERVAL_LENGTH, 0);  

    float prevLogEstimate = 0;
    uint32_t nextLogEstimate = 1;

    while (smooth.size() < factorBase.size() + 20) {
       
   
        for (uint32_t i = 1, a = intervalStart + 1; i < INTERVAL_LENGTH; ++i, ++a) {
            if (nextLogEstimate <= a) {
                const mpz_class Q = (a + sqrtN) * (a + sqrtN) - N;
                prevLogEstimate = mpz_sizeinbase(Q.get_mpz_t(), 2);    // ~log_2(Q)
                nextLogEstimate = nextLogEstimate * 1.8 + 1;
            }
            logApprox[i] = prevLogEstimate;
        }

      
      
        for (uint32_t i = 0; i < factorBase.size(); ++i) {
            const uint32_t p = factorBase[i];
            const float logp = std::log(factorBase[i]) / std::log(2);

            // Sieve first sequence.
            while (startIndex.first[i] < intervalEnd) {
                logApprox[startIndex.first[i] - intervalStart] -= logp;
                startIndex.first[i] += p;
            }

            if (p == 2)
                continue; 

            // Sieve second sequence.
            while (startIndex.second[i] < intervalEnd) {
                logApprox[startIndex.second[i] - intervalStart] -= logp;
                startIndex.second[i] += p;
            }
        }
     
        const float threshold = std::log(factorBase.back()) / std::log(2);
        for (uint32_t i = 0, a = intervalStart; i < INTERVAL_LENGTH; ++i, ++a) {
            if (std::fabs(logApprox[i]) < threshold) {
                mpz_class Q = (a + sqrtN) * (a + sqrtN) - N;
                std::vector<uint32_t> factors;

                for (uint32_t j = 0; j < factorBase.size(); ++j) {
                    const uint32_t p = factorBase[j];
                    while (mpz_divisible_ui_p(Q.get_mpz_t(), p)) {
                        mpz_divexact_ui(Q.get_mpz_t(), Q.get_mpz_t(), p);
                        factors.push_back(j);                    }
                }
                if (Q == 1) {
                    smoothFactors.push_back(factors);
                    smooth.push_back(a);
                }

                if (smooth.size() >= factorBase.size() + 20)
                    break;
            }
        }
      
        intervalStart += INTERVAL_LENGTH;
        intervalEnd += INTERVAL_LENGTH;
    }


 

    /*
     * Step 3.1
     *
     * Construct a binary matrix M with M_ij = the parity of the i:th prime factor
     * from the factor base in the factorization of the j:th B-smooth number.
     */
    Matrix M(factorBase.size(), smoothFactors.size() + 1);
    for (uint32_t i = 0; i < smoothFactors.size(); ++i) {
        for (uint32_t j = 0; j < smoothFactors[i].size(); ++j) {
            M(smoothFactors[i][j], i).flip();
        }
    }

    /*
     * Step 3.2
     *
     * Reduce the matrix to row echelon form and solve it repeatedly until a factor
     * is found.
     */
    M.reduce();
    mpz_class a;
    mpz_class b;

    do {
        std::vector<uint32_t> x = M.solve();

        a = 1;
        b = 1;

        /*
         * Calculate b = product(smooth[i] + sqrt(N)).
         *
         * Also calculate the the power of each prime in a's decomposition on the
         * factor base.
         */
        std::vector<uint32_t> decomp(factorBase.size(), 0);
        for (uint32_t i = 0; i < smoothFactors.size(); ++i) {
            if (x[i] == 1) {
                for(uint32_t p = 0; p < smoothFactors[i].size(); ++p)
                    ++decomp[smoothFactors[i][p]];
                b *= (smooth[i] + sqrtN);
            }
        }

        for(uint32_t p = 0; p < factorBase.size(); ++p) {
            for(uint32_t i = 0; i < (decomp[p] / 2); ++i)
                a *= factorBase[p];
        }

    } while (a % N == b % N || a % N == (- b) % N + N);

    mpz_class factor;
    mpz_gcd(factor.get_mpz_t(), mpz_class(b - a).get_mpz_t(), N.get_mpz_t());

    return factor;
}
