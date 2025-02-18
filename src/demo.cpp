#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-remez.h"
#include "ranking.h"
#include "ranking_compos.h"
#include "optimized_degrees.h"
#include "minimum.h"
#include "sorting.h"
#include "openfhe.h"
#include "pke/cryptocontext.h"
#include "key/privatekey.h"
#include "key/publickey.h"
#include "math/chebyshev.h"
#include "schemerns/rns-scheme.h"
#include "scheme/ckksrns/ckksrns-cryptoparameters.h"
#include <iostream>
#include <NTL/RR.h>
#include <cmath>
#include "optimized_degrees.h"
#include <vector>
#include <functional>
#include <cstdlib>  // For system()
#include <string>   // For std::string and std::to_string
#include <sstream>  // For std::ostringstream
#include <utility>  // For std::pair
#include <fstream>    // For std::ifstream
#include <filesystem>
#include <iostream>
#include <unistd.h>
#include <cstdio>
#include <memory>
#include <string>

using namespace std;
using namespace NTL;

#include <cmath>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iomanip>

// Placeholder functions for MP() and ME()
double MP(const std::pair<double, double>& R_squigly, int M_deg) {
    // Implement the MP function here
    // For now, it returns a dummy value
    return 0.0;
}
double ME(const std::pair<double, double>& R_squigly, int M_deg) {
    // Implement the ME function here
    // For now, it returns a dummy value
    return 0.0;
}
// Function to evaluate a polynomial at a given point
double evaluate_polynomial(const std::vector<double>& coefficients, double x) {
    double result = 0.0;
    double power_of_x = 1.0;
    for (double coeff : coefficients) {
        result += coeff * power_of_x;
        power_of_x *= x;
    }
    return result;
}

template <typename T>
int sign(T value) {
    if (value > 0) {
        return 1;
    } else if (value < 0) {
        return -1;
    } else {
        return 0;
    }
}

template <typename T>
int abs(T value) {
    if (value < 0) {
        return -(value);
    } else {
        return value;
    }
}

std::vector<int> calculateRanks(const std::vector<double>& vec) {
    int n = vec.size();
    
    // Create a vector of indices from 0 to vec.size()-1
    std::vector<int> indices(n);
    for (int i = 0; i < n; ++i) {
        indices[i] = i;
    }
    
    // Sort the indices based on the corresponding values in vec
    std::sort(indices.begin(), indices.end(), [&vec](int a, int b) {
        return vec[a] < vec[b];  // Sort indices by comparing the values in vec
    });
    
    // Create a rank vector
    std::vector<int> ranks(n);
    
    // Assign ranks based on sorted indices
    for (int i = 0; i < n; ++i) {
        ranks[indices[i]] = i + 1;  // Rank starts from 1
    }
    
    return ranks;
}


void checkMemoryUsage(pid_t pid) {
    std::string command = "ps -o pid,rss,cmd -p " + std::to_string(pid);
    char buffer[128];
    std::string result;

    // Open the command for reading
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
        std::cerr << "Failed to run ps command." << std::endl;
        return;
    }

    // Read the output
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result += buffer;
    }

    pclose(pipe);

    // Print the output
    std::cout << result << std::endl;
}

namespace fs = std::filesystem;

std::string getCurrentDate() {
    std::time_t t = std::time(nullptr);
    std::tm tm{};
#ifdef _WIN32
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char buffer[11];
    std::strftime(buffer, sizeof(buffer), "%Y%m%d", &tm);
    return buffer;
}

void clearFileContents(const fs::path& filePath) {
    std::ofstream ofs(filePath, std::ios::trunc); // Open file in truncate mode
    if (!ofs) {
        std::cerr << "Error clearing file: " << filePath << std::endl;
    }
}

void backup_text() {
    fs::path sourceDir = "../text";
    std::cout << "Backing up files in 'text' directory...\n";

    // Check if "text" directory exists
    if (!fs::exists(sourceDir) || !fs::is_directory(sourceDir)) {
        std::cout << "Folder 'text' does not exist.\n";
        return;
    }

    // Check if there are files inside the directory
    bool hasFiles = false;
    for (const auto& entry : fs::directory_iterator(sourceDir)) {
        if (fs::is_regular_file(entry)) {
            hasFiles = true;
            break;
        }
    }

    if (!hasFiles) {
        std::cout << "No files found in 'text'.\n";
        return;
    }

    // Create a new backup folder with today's date
    std::string dateStr = getCurrentDate();
    fs::path backupDir = "text" + dateStr;

    if (!fs::exists(backupDir)) {
        fs::create_directory(backupDir);
    }

    // Copy files and subdirectories to the new backup folder
    for (const auto& entry : fs::recursive_directory_iterator(sourceDir)) {
        fs::path destPath = backupDir / entry.path().lexically_relative(sourceDir);

        try {
            if (fs::is_directory(entry)) {
                fs::create_directory(destPath);
            } else if (fs::is_regular_file(entry)) {
                fs::copy(entry.path(), destPath, fs::copy_options::overwrite_existing);
            }
        } catch (const fs::filesystem_error& e) {
            std::cerr << "Error copying: " << e.what() << std::endl;
        }
    }

    // Clear contents of all files in "text"
    for (const auto& entry : fs::directory_iterator(sourceDir)) {
        if (fs::is_regular_file(entry)) {
            clearFileContents(entry.path());
        }
    }

    std::cout << "Backup completed to '" << backupDir.string() << "' and original files cleared.\n";
}

int main()
{

    std::cout << std::endl <<
        "Test: computing time consumption and error for Minimax composite Approximation"
        << std::endl << std::endl;

    ////////////////////////////////////////////////////////////////////////
    //                       Setting up variables                         //
    ////////////////////////////////////////////////////////////////////////

    std::cout << "Setting up variables... \n" << std::endl;
    
    double lowerBounddepth = -1;
    double upperBounddepth = 1;
    float alpha = 6.64385619;
    std::string alpha_string = "6.64385619";
    
    // Do we want to do new measurements fo CC?
    bool doCC = false;

    std::string output_directory = "../remez_outputs_cluster";
    
    std::cout << std::fixed << std::setprecision(9);
    std::cout << "Alpha: " << alpha << std::endl;
    std::cout << "E: " << alpha2error(alpha) << std::endl;

    std::vector<int> degrees;

    // generate vector
    int num_points = 8;
    double start = alpha2error(alpha)/2;
    double end = 1.0;
    std::vector<double> v(num_points);

    // Calculate step size
    double step = (end - start) / (num_points - 1);

    // Lambda to generate each point
    double current = start; // variable to hold the current point
    std::generate_n(v.begin(), num_points, [&]() mutable {
        double value = current;
        current += step;
        return value;
    });

    std::cout << "---------------------------------\n" << std::endl;
    
    ////////////////////////////////////////////////////////////////////////
    //                       Setting up CKKS scheme                       //
    ////////////////////////////////////////////////////////////////////////

    std::cout << "Setting up CKKS scheme \n" << std::endl;


    const size_t vectorLength = 8;
    const usint integralPrecision       = 12;
    const usint decimalPrecision        = 44;
    const usint multiplicativeDepth     = alpha2depthCompos(alpha) + 15;
    const usint numSlots                = vectorLength * vectorLength;
    const bool enableBootstrap          = false;
    
    const usint ringDim                 = 0;
    const bool verbose                  = true;

    CryptoContext<DCRTPoly> cryptoContext = generateCryptoContext(
            integralPrecision,
            decimalPrecision,
            multiplicativeDepth,
            numSlots,
            enableBootstrap,
            ringDim,
            verbose
    );

    // Generating public/private key pair, relinearization, and rotation keys
    std::vector<int32_t> indices = getRotationIndices(vectorLength);
    KeyPair<DCRTPoly> keyPair = keyGeneration(
            cryptoContext,
            indices,
            numSlots,
            enableBootstrap,
            verbose
    );

    cryptoContext->EvalMultKeyGen(keyPair.secretKey);

    Ciphertext<DCRTPoly> vC = cryptoContext->Encrypt(
        keyPair.publicKey,
        cryptoContext->MakeCKKSPackedPlaintext(v)
    );
    std::cout << "---------------------------------\n" << std::endl;
    
      //////////////////////////////////////////////////////////////////////////////////////
     //                            CC measurements for degrees                           //
    //////////////////////////////////////////////////////////////////////////////////////
    
    if (!doCC) {
        std::cout << "Skipping CC measurements \n" << std::endl;
    }

    if (doCC) {
        std::cout << "Measuring evaluation times for calculating optimal degrees" << std::endl;
        ConstCiphertext<DCRTPoly> constCiphertext = vC;
        double lowerBound = 0;
        double upperBound = 10;

        // Function we want to approximate
        std::function<double(double)> func = [](double x) {
        return (x > 0) - (x < 0);
        };

        
        
        // Calculate evaluation times for polynomials of different degrees
        for(int32_t polyDegree=3; polyDegree<=63; polyDegree+=2){ 
            
            std::cout << "working on degree: " << polyDegree << std::endl;

            // Calculate coefficients from Chebyshev approximation
            std::vector<double> coefficients = EvalChebyshevCoefficients(func, lowerBound, upperBound, polyDegree);

            int n_max = 30;
            std::vector<int> timeVector;


            int logBase2Degree = std::ceil(log2(polyDegree));
            
            // Loop over ciphertext of different depths
            for(int i=n_max; i>=0; i--){
                auto ciphertextReduced = vC;
                cryptoContext->EvalMultKeyGen(keyPair.secretKey);

                for(int j=n_max; j>=0; j--)
                {   
                    if(j > i)
                    { 
                        timeVector.insert(timeVector.begin(), 10000000);
                    }
                    else if(j < logBase2Degree) {
                        timeVector.insert(timeVector.begin(), 10000000);
                    }
                    else 
                    {
                        auto t1 = std::chrono::high_resolution_clock::now();
                        auto result = cryptoContext->EvalChebyshevSeries(ciphertextReduced, coefficients, lowerBound, upperBound);
                        auto t2 = std::chrono::high_resolution_clock::now();
                        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1);
                        timeVector.insert(timeVector.begin(), ms_int.count());
                        // timeVector.push_back(ms_int);
                    }
                    if (j == 0){
                        break;
                    }
                    ciphertextReduced = cryptoContext->EvalMult(ciphertextReduced, 1);
                }
            }


            ////////////////////////////////////////////////////////////////////////
            //                         Write measurements to file                 //
            ////////////////////////////////////////////////////////////////////////
            

            // if the ../text/ folder is already filled with the results, 
            // create a copy of the text file and save it in a folder named "textYYYYMMDD" (using the current date) folder
            backup_text();

            // Set output file path
            std::ostringstream oss;
            oss << "../text/CC" << polyDegree << ".txt";
            std::string path = oss.str();

            std::ofstream outFile(path);

            // Check if the file opened successfully
            if (!outFile) {
            std::cerr << "Error opening file for writing" << std::endl;
            return 1;
            }
            std::cout << "writing results: " << std::endl;
            // Iterate through the vector and write each value to the file
            for (const int& value : timeVector) {
                outFile << value << "\n";  // Write each value followed by a newline
            }
            std::cout << "done writing results \n" << std::endl;

            outFile.close();
        }
    }

    std::cout << "---------------------------------\n" << std::endl;

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     //      approximate using multi-interval remez and optimal degrees and calculate the optimal degrees for alpha      //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Construct the path
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(9);
    oss << output_directory << "/" << alpha_string;
    std::string coefficient_path = oss.str();
    
    // Extract the directory part of the path
    std::filesystem::path path(coefficient_path);
    std::filesystem::path directory = path.parent_path();
    std::cout << "Checking if the directory exists: " << path << std::endl;

    degrees = alpha2degreeCompos(alpha);
    std::cout << "Degrees: " << degrees << std::endl;

    // check if degrees are part of precalculated polynomials
    if (size(degrees) > 0) {
        std::cout << "Using known degrees \n" << std::endl;
    } else {
        std::cout << "Calculating optimal degrees for alpha..." << std::endl;
        std::cout << "---------------------------------\n" << std::endl;

        long maxdeg = 63;		// 31 or 63
        bool is_comp = true;	// true: comparison operation, false: max/ReLU operation
        long max_factor = 1;	// max_factor = 1 for comparison operation. max_factor > 1 for max/ReLU operation 
        RR epsilon = max_factor * pow(RR(2),RR(-alpha));
        
        
        std::vector<int> degrees = compute_min_multdepth(RR(alpha), epsilon, maxdeg, is_comp);
        
        std::cout << "Degrees: " << degrees << std::endl;
        std::cout << "---------------------------------\n" << std::endl;
    }
    

    // Check if the directory exists, no need to recalculate coefficients
    if (std::filesystem::exists(path)) {
        // check if the output directory exists 
        std::cout << "Using previous approximations \n" << std::endl;
    } else {
        std::cout << "Approximation" << std::endl;
        std::cout << "---------------------------------\n" << std::endl;

        std::cout << "Creating new approximations " << std::endl;
        // Construct the command string
        // Convert the vector to a space-separated string
        std::ostringstream oss_degrees;
        for (size_t i = 0; i < degrees.size(); ++i) {
            oss_degrees << degrees[i];
            if (i < degrees.size() - 1) {
                oss_degrees << " ";
            }
        }
        std::string degrees_str = oss_degrees.str();
        std::string command = "sudo python3 ../remez_approximation_pipeline.py " + std::to_string(alpha) + " " + degrees_str;
        
        // Print the command for debugging
        std::cout << "Running command: " << command << std::endl;
        
        // Run the Python script
        int result = std::system(command.c_str());
        
        if (result == 0) {
            std::cout << "Approximation script executed successfully!" << std::endl;
        } else {
            std::cerr << "Error: Approximation script execution failed with code " << result << std::endl;
        }

        std::cout << "---------------------------------\n" << std::endl;
    }

    
      //////////////////////////////////////////////////////////////////////////////////////
     //          Computing Time Consumption Composite Minimax Polynomial                 //
    //////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Computing Time Consumption Composite Minimax Polynomial..." << std::endl;

    // Do measurements
    int ms_int = 0;
    auto resultC = vC;
    std::vector<std::vector<double>> coefficient_list;
    for(size_t i = 0; i < degrees.size(); i++) {
        
        // Setup the path to the coefficients file
        auto polyDegree = degrees[i];
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(9);
        oss << output_directory << "/" << alpha_string << "/coefficients_" << polyDegree << "_" << i+1 << ".txt";
        std::string coefficient_path = oss.str();
        
        // Read float values from the coefficients file
        std::ifstream infile(coefficient_path);
        std::vector<float> coefficients;
        float value;
        while (infile >> value) {
            coefficients.push_back(value);
        }
        infile.close();

        // Convert float values to double
        std::vector<double> double_coefficients;
        double_coefficients.reserve(coefficients.size()); // Reserve space for efficiency
        for (float value : coefficients) {
            double_coefficients.push_back(static_cast<double>(value));
        }
        coefficient_list.push_back(double_coefficients);
        
        // Run the sign function approximation over the encrypted vector
        auto t1 = std::chrono::high_resolution_clock::now();
        resultC = cryptoContext->EvalChebyshevSeries(resultC, double_coefficients, lowerBounddepth, upperBounddepth);
        auto t2 = std::chrono::high_resolution_clock::now();
        ms_int += std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        
        std::reverse(double_coefficients.begin(), double_coefficients.end());\

    }
    
    // print results
    std::cout << "Time taken: " << ms_int << std::endl;
    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Results            : " << result << "\n" << std::endl;
    
    // Calculate max and average error 
    float max_error = 0.0;
    float error = 0.0;
    for (size_t i = 0; i < result.size(); ++i) {
        if (max_error < float(abs(result[i] - sign(v[i])))) {
            max_error = float(abs(result[i] - sign(v[i])));
        }
        error += float(abs(result[i] - sign(v[i])));
    }
    std::cout.precision(14);
    std::cout << "Average error     : " << (float(error / float(v.size()))) << std::endl;
    std::cout << "Max error         : " << max_error  << "\n" << std::endl;
    
      //////////////////////////////////////////////////////////////////////////////////////////////
     //    Computing Time Consumption of Ranking Using Composite Minimax Polynomial              //
    //////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Computing time consumption of ranking algorithm using Composite Minimax Polynomial..." << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    // Run ranking algorithm using 
    resultC = rankComposCipher(
        vC,
        coefficient_list,
        vectorLength,
        -1.0, 1.0,
        degrees,
        alpha,
        alpha_string
    );
    auto t2 = std::chrono::high_resolution_clock::now();
    ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();

    std::cout << "Time taken: " << ms_int << std::endl;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    result = resultP->GetRealPackedValue();
    std::cout.precision(3);
    std::cout << "Ranks            : " << result << "\n" << std::endl;

    std::vector<int> ranks = calculateRanks(v);

    error = 0.0;
    for (size_t i = 0; i < result.size(); ++i) {
        error += float(abs(result[i] - ranks[i]));
    }
    
    std::cout << "Average error     : " << (float(error / float(v.size()))) << "\n" << std::endl;
    
    
      //////////////////////////////////////////////////////////////////////////////////////
     //           Computing Time Consumption Minimum Using Min function                  //
    //////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "Computing Time Consumption for Calculating Minimum Using The Min Function..." << std::endl;

    const usint indicatorDepth = 11;
    int indic_degree = depth2degree(indicatorDepth);
    double min_dif_bound = (v[1] - v[0]);


    resultC = min_adapted_max(
             vC,
             vectorLength,
             -1.0, 1.0,
             245,
             indic_degree,
             coefficient_list,
             degrees,
             min_dif_bound
    );

    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    result = resultP->GetRealPackedValue();
    
    double computedMin = 0.0;
    double norm = 0.0;
    for (size_t i = 0; i < v.size(); i++)
    {
        computedMin += v[i] * result[i];
        norm += result[i];
    }
    computedMin /= norm;
    std::cout << "Minimum           : " << computedMin << std::endl;
    std::cout << "Expected minimum  : " << min(v) << std::endl << std::endl;

    return 0;

}
