#include "ranking_compos.h"

using namespace std;
using namespace NTL;


Ciphertext<DCRTPoly> rankComposCipher(
    Ciphertext<DCRTPoly> c,
    std::vector<std::vector<double>> coefficient_list,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    std::vector<int> M_degs,
    float alpha,
    std::string alpha_string,
    bool cmpGt
)
{
    if (!cmpGt)
    {
        c = compare_composite(
            replicateRow(c, vectorLength),
            replicateColumn(transposeRow(c, vectorLength, true), vectorLength),
            coefficient_list,
            leftBoundC, rightBoundC, M_degs, alpha, alpha_string
        );
    }
    else
    {
        c = compareGt(
            replicateRow(c, vectorLength),
            replicateColumn(transposeRow(c, vectorLength, true), vectorLength),
            leftBoundC, rightBoundC, M_degs[0],
            0.01
        );
    }

    c = sumRows(c, vectorLength);
    c = c + (!cmpGt ? 0.5 : 1.0);

    return c;
}



std::vector<Ciphertext<DCRTPoly>> rankComposCipherVector(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    bool cmpGt,
    bool complOpt
)
{
    const size_t numCiphertext = c.size();

    static std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;

    std::cout << "===================================\n";
    std::cout << "Replicate\n";
    std::cout << "===================================\n";

    std::vector<Ciphertext<DCRTPoly>> replR(numCiphertext);
    std::vector<Ciphertext<DCRTPoly>> replC(numCiphertext);

    start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for collapse(2)
    for (size_t loopID = 0; loopID < 2; loopID++)
    {
        for (size_t j = 0; j < numCiphertext; j++)
        {
            if (loopID == 0)
            {
                #pragma omp critical
                {std::cout << "ReplicateRow - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                replR[j] = replicateRow(c[j], subVectorLength);

                #pragma omp critical
                {std::cout << "ReplicateRow - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
            else
            {
                #pragma omp critical
                {std::cout << "ReplicateColumn - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                replC[j] = replicateColumn(transposeRow(c[j], subVectorLength, true), subVectorLength);

                #pragma omp critical
                {std::cout << "ReplicateColumn - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
        }
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
    if (!complOpt)
    {
        std::cout << "===================================\n";
        std::cout << "Compare\n";
        std::cout << "===================================\n";

        std::vector<Ciphertext<DCRTPoly>> C(numCiphertext);
        std::vector<bool> Cinitialized(numCiphertext, false);

        start = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for collapse(2)
        for (size_t j = 0; j < numCiphertext; j++)
        {
            for (size_t k = 0; k < numCiphertext; k++)
            {
                #pragma omp critical
                {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                Ciphertext<DCRTPoly> Cjk;
                if (!cmpGt)
                {
                    Cjk = compare(
                        replR[j],
                        replC[k],
                        leftBoundC, rightBoundC, degreeC
                    );
                }
                else
                {
                    Cjk = compareGt(
                        replR[j],
                        replC[k],
                        leftBoundC, rightBoundC, degreeC,
                        0.01
                    );
                }

                #pragma omp critical
                {
                if (!Cinitialized[j]) { C[j] = Cjk; Cinitialized[j] = true; }
                else                  { C[j] = C[j] + Cjk;                  }
                }
                
                #pragma omp critical
                {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed_seconds = end - start;
        std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
                elapsed_seconds.count() << "s)" << std::endl << std::endl;

        std::cout << "===================================\n";
        std::cout << "Sum\n";
        std::cout << "===================================\n";

        std::vector<Ciphertext<DCRTPoly>> s(numCiphertext);

        start = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for
        for (size_t j = 0; j < numCiphertext; j++)
        {
            #pragma omp critical
            {std::cout << "SumRows - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

            s[j] = sumRows(C[j], subVectorLength) + (!cmpGt ? 0.5 : 1.0);
            
            #pragma omp critical
            {std::cout << "SumRows - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed_seconds = end - start;
        std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
                elapsed_seconds.count() << "s)" << std::endl << std::endl;
        
        return s;
    }
    else
    {
        std::cout << "===================================\n";
        std::cout << "Compare\n";
        std::cout << "===================================\n";

        std::vector<Ciphertext<DCRTPoly>> Cv(numCiphertext);
        std::vector<Ciphertext<DCRTPoly>> Ch(numCiphertext);
        std::vector<bool> Cvinitialized(numCiphertext, false);
        std::vector<bool> Chinitialized(numCiphertext, false);

        start = std::chrono::high_resolution_clock::now();

        const size_t numReqThreads = numCiphertext * (numCiphertext + 1) / 2;
        std::cout << "Number of required threads: " << numReqThreads << std::endl;

        // for (size_t j = 0; j < numCiphertext; j++)
        // {
        //     for (size_t k = j; k < numCiphertext; k++)
        //     {
        // Collapse(2) with two nested for-loops creates issues here.
        #pragma omp parallel for
        for (size_t i = 0; i < numReqThreads; i++)
        {
            // Computing the indeces
            size_t j, k, counter = 0;
            k = 0;
            bool loopCond = true;
            for (j = 0; j < numCiphertext && loopCond; j++)
                for (k = j; k < numCiphertext && loopCond; k++)
                    if (counter++ == i) loopCond = false;
            j--; k--;

            #pragma omp critical
            {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

            Ciphertext<DCRTPoly> Cjk;
            if (!cmpGt)
            {
                Cjk = compare(
                    replR[j],
                    replC[k],
                    leftBoundC, rightBoundC, degreeC
                );
            }
            else
            {
                Cjk = compareGt(
                    replR[j],
                    replC[k],
                    leftBoundC, rightBoundC, degreeC,
                    0.001
                );
            }

            #pragma omp critical
            {
            if (!Cvinitialized[j]) { Cv[j] = Cjk; Cvinitialized[j] = true; }
            else                   { Cv[j] = Cv[j] + Cjk;                  }
            }

            if (j != k)
            {
                Ciphertext<DCRTPoly> Ckj = 1.0 - Cjk;

                #pragma omp critical
                {
                if (!Chinitialized[k]) { Ch[k] = Ckj; Chinitialized[k] = true; }
                else                   { Ch[k] = Ch[k] + Ckj;                  }
                }
            }
            
            #pragma omp critical
            {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed_seconds = end - start;
        std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
                elapsed_seconds.count() << "s)" << std::endl << std::endl;

        std::cout << "===================================\n";
        std::cout << "Sum\n";
        std::cout << "===================================\n";

        std::vector<Ciphertext<DCRTPoly>> sv(numCiphertext);
        std::vector<Ciphertext<DCRTPoly>> sh(numCiphertext);
        std::vector<Ciphertext<DCRTPoly>> s(numCiphertext);
        std::vector<bool> sinitialized(numCiphertext, false);

        start = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for collapse(2)
        for (size_t loopID = 0; loopID < 2; loopID++)
        {
            for (size_t j = 0; j < numCiphertext; j++)
            {
                if (loopID == 0)
                {
                    #pragma omp critical
                    {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                    sv[j] = sumRows(Cv[j], subVectorLength, true);
                    
                    #pragma omp critical
                    {
                    if (!sinitialized[j]) { s[j] = sv[j] + (!cmpGt ? 0.5 : 1.0); sinitialized[j] = true; }
                    else                  { s[j] = s[j] + sv[j];                                             }
                    }

                    #pragma omp critical
                    {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
                }
                else
                {
                    #pragma omp critical
                    {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                    if (j > 0)
                    {
                        sh[j] = sumColumns(Ch[j], subVectorLength, true);
                        sh[j] = transposeColumn(sh[j], subVectorLength);

                        #pragma omp critical
                        {
                        if (!sinitialized[j]) { s[j] = sh[j] + (!cmpGt ? 0.5 : 1.0); sinitialized[j] = true; }
                        else                  { s[j] = s[j] + sh[j];                                             }
                        }
                    }

                    #pragma omp critical
                    {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
                }
            }
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed_seconds = end - start;
        std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
                elapsed_seconds.count() << "s)" << std::endl << std::endl;
        
        return s;
    }
    
}


// fractional rank
std::vector<double> rankComposPlain(
    const std::vector<double> &vec,
    double epsilon
)
{
    std::vector<double> ranking(vec.size(), 0.5);
    for (size_t i = 0; i < vec.size(); i++)
        for (size_t j = 0; j < vec.size(); j++)
        { 
            if (std::abs(vec[i] - vec[j]) <= epsilon)
                ranking[i] += 0.5;
            else if (vec[i] > vec[j])
                ranking[i] += 1.0;
        }

    return ranking;
}
