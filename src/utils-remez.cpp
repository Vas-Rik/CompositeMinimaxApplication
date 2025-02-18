#include "utils-remez.h"

void run_remez(float approximation_param, std::vector<std::pair<RR, RR>> D, std::vector<int> degrees) {
     
    // Convert float and int to strings
    std::string approximation_param_str = std::to_string(approximation_param);
    // std::string degrees_str = std::to_string(degrees);
    
    // Convert the list of degrees to a string
    std::ostringstream Degree_stream;
    Degree_stream << "[";
    for (const auto& deg : degrees) {
        Degree_stream << deg << ",";
        // std::cout << deg << ",";
    }
    Degree_stream << "] ";
    // std::cout << endl


    std::string degrees_str = Degree_stream.str();
    degrees_str.pop_back();  // Remove the trailing space

    // Convert the list of tuples to a string
    std::ostringstream D_stream;
    D_stream << "[";
    for (const auto& pair : D) {
        D_stream << "(" << pair.first << "," << pair.second << "),";
        // std::cout << "(" << pair.first << "," << pair.second << ")\n";
    }
    D_stream << "] ";
    std::string D_str = D_stream.str();
    D_str.pop_back();  // Remove the trailing space

    // Construct the command string
    std::string command = "python3 ../src/remez_2.py " + approximation_param_str + " \"" + D_str + "\" " + degrees_str;

    std::cout << command << std::endl;

    // Execute the command using system()
    int result = system(command.c_str());

    // Check the result of the command execution
    if (result == 0) {
        std::cout << "Command executed successfully.\n \n";
    } else {
        std::cerr << "Command failed.\n";
    }
}



std::vector<double> read_coefficients(int degree, int i) {
    // File name to read
    std::ostringstream oss;
    oss << "../remez_outputs/coefficients_" << degree << "_" << i << ".txt";
    std::string filename = oss.str();

    // Open the file
    std::ifstream infile(filename);

    // Check if the file was opened successfully
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << "\n";
    }

    // Vector to store the float values
    std::vector<double> coeffs;

    // Read the file line by line
    std::string line;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        double value;
        
        // Attempt to convert the line to a float
        if (ss >> value) {
            coeffs.push_back(value);  // Store the float in the vector
        } else {
            std::cerr << "Error: Invalid float value on line: " << line << "\n";
        }
    }

    // Close the file
    infile.close();

    return coeffs;
}


std::vector<std::vector<double>> read_coefficients_alpha(float alpha, std::vector<int> degrees) {
    
    std::vector<std::vector<double>> coefficient_list = {};

    int i = 1;
    for (const int& degree : degrees) {
        std::ostringstream oss;
        oss << "../coefficients_minimax/" << alpha << "/" << degree << "_" << i << ".txt";
        std::string filename = oss.str();

        // Open the file
        std::ifstream infile(filename);

        // Check if the file was opened successfully
        if (!infile.is_open()) {
            std::cerr << "Error: Could not open the file " << filename << "\n";
        }

        // Vector to store the float values
        std::vector<double> coeffs;

          // Read the file line by line
        std::string line;
        while (std::getline(infile, line)) {
            std::stringstream ss(line);
            double value;
            
            // Attempt to convert the line to a float
            if (ss >> value) {
                coeffs.push_back(value);  // Store the float in the vector
            } else {
                std::cerr << "Error: Invalid float value on line: " << line << "\n";
            }
        }

        // Close the file
        infile.close();
        coefficient_list.push_back(coeffs);
        i+=1;
    }

    return coefficient_list;
}



float read_error() {
    // File name to read
    std::string filename = "../remez_outputs/error.txt";

    // Open the file
    std::ifstream infile(filename);

    // Check if the file was opened successfully
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << "\n";
    }

    // Vector to store the float values
    float error = 0.0;

    // Read the file line by line
    std::string line;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        float value;
        // Attempt to convert the line to a float
        if (ss >> value) {
            error = value;  // Store the float in the vector
        } else {
            std::cerr << "Error: Invalid float value on line: " << line << "\n";
        }
    }

    // Close the file
    infile.close();

    // Output the values to verify (optional)
    std::cout << "\n Stored error value :\n";
    std::cout << error << "\n";
    
    return error;
}


// Function to compute the Chebyshev polynomial of the first kind T_n(x)
float chebyshevT(int n, float x) {
    if (n == 0) return 1.0f;
    if (n == 1) return x;
    float T0 = 1.0f;
    float T1 = x;
    float Tn = 0;
    for (int i = 2; i <= n; ++i) {
        Tn = 2.0f * x * T1 - T0;
        T0 = T1;
        T1 = Tn;
    }
    return Tn;
}


// Function to evaluate the polynomial at point x using Chebyshev basis
float evaluateChebyshevBasis2(const std::vector<double>& coefficients, float x) {
    float result = 0.0f;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        result += static_cast<float>(coefficients[i]) * chebyshevT(static_cast<int>(i), x);
    }
    return result;
}

float sign(float x) {
    if      (x > 0)       return 1;
    if      (x < 0)       return -1;
    else                  return 0;
}

// Function to evaluate a linear combination of Chebyshev polynomials of the first kind
float evaluateChebyshevBasis(const std::vector<double>& coefficients, float x) {
    int n = coefficients.size();
    double doub_x = static_cast<double>(x);
    double c0 = 0;
    double c1 = 0; 
    if (n == 1) {
        c0 = coefficients[0];
        c1 = 0;
    } else if (n == 2) {
        c0 = coefficients[0];
        c1 = coefficients[1];
    } else {
        double x2 = 2 * doub_x;
        c0 = coefficients[n-2];
        c1 = coefficients[n-1];
        for (int i = 3; i <= n; i++) {
            double tmp = c0;
            c0 = coefficients[-i] - c1;
            c1 = tmp + c1 * x2;
        }
    }
    return static_cast<float>(c0 + c1 * doub_x);
}


// Function to check if max error is less than err
float max_error_interval(const std::vector<double>& coefficients, float err) {
    float max_error = 0.0;
    const int num_samples = 1000;  // Number of samples to evaluate polynomial
    float step = (1 - (err / 2)) / num_samples;

    // Sample points in the interval [err/2, 1]
    for (int i = 0; i <= num_samples; ++i) {
        float x = (err / 2) + (i * step);
        float poly_value = evaluateChebyshevBasis(coefficients, x);
        max_error = std::max(max_error, std::abs(poly_value - sign(x)));
    }

    // Return true if max error is less than the given err
    return max_error;
}

