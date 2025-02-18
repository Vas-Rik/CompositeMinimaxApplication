#include "utils-eval.h"



usint depth2degree(
    const usint depth
)
{
    switch(depth)
    {
        case 4:     return 5;
        case 5:     return 13;
        case 6:     return 27;
        case 7:     return 59;
        case 8:     return 119;
        case 9:     return 247;
        case 10:    return 495;
        case 11:    return 1007;
        case 12:    return 2031;

        case 13:    return 4031;
        case 14:    return 8127;
        default:    return -1;
    }
}



usint alpha2degreeOpenFHE(const float alpha)
{
    const float epsilon = 0.0001f; // Tolerance for comparing floating-point values
    
    if (std::fabs(alpha - 3.091503713) < epsilon)       return 3;
    else if (std::fabs(alpha - 3.62961696) < epsilon)   return 5;
    else if (std::fabs(alpha - 4.026158059) < epsilon)  return 7;
    else if (std::fabs(alpha - 4.339009337) < epsilon)  return 9;
    else if (std::fabs(alpha - 4.596911194) < epsilon)  return 11;
    else if (std::fabs(alpha - 4.816362316) < epsilon)  return 13;
    else if (std::fabs(alpha - 5.006872028) < epsilon)  return 15;
    else if (std::fabs(alpha - 5.175246808) < epsilon)  return 17;
    else if (std::fabs(alpha - 5.326378453) < epsilon)  return 19;
    else if (std::fabs(alpha - 5.449221655) < epsilon)  return 21;
    else if (std::fabs(alpha - 5.575712623) < epsilon)  return 23;
    else if (std::fabs(alpha - 5.69235669) < epsilon)   return 25;
    else if (std::fabs(alpha - 5.800459805) < epsilon)  return 27;
    else if (std::fabs(alpha - 5.900865483) < epsilon)  return 29;
    else if (std::fabs(alpha - 5.99483256) < epsilon)   return 31;
    else if (std::fabs(alpha - 6.082714136) < epsilon)  return 33;
    else if (std::fabs(alpha - 6.166037129) < epsilon)  return 35;
    else if (std::fabs(alpha - 6.244093359) < epsilon)  return 37;
    else if (std::fabs(alpha - 6.318384415) < epsilon)  return 39;
    else if (std::fabs(alpha - 6.389465946) < epsilon)  return 41;
    else if (std::fabs(alpha - 6.448095937) < epsilon)  return 43;
    else if (std::fabs(alpha - 6.510721532) < epsilon)  return 45;
    else if (std::fabs(alpha - 6.573560838) < epsilon)  return 47;
    else if (std::fabs(alpha - 6.63309892) < epsilon)   return 49;
    else if (std::fabs(alpha - 6.690314796) < epsilon)  return 51;
    else if (std::fabs(alpha - 6.745604685) < epsilon)  return 53;
    else if (std::fabs(alpha - 6.799211286) < epsilon)  return 55;
    else if (std::fabs(alpha - 6.849537749) < epsilon)  return 57;
    else if (std::fabs(alpha - 6.899942561) < epsilon)  return 59;
    else if (std::fabs(alpha - 6.946884719) < epsilon)  return 61;
    else if (std::fabs(alpha - 6.993709998) < epsilon)  return 63;
    else if (std::fabs(alpha - 8.939735521) < epsilon)  return 247;
    else if (std::fabs(alpha - 9.939149282) < epsilon)  return 495;
    else if (std::fabs(alpha - 10.962565805) < epsilon)  return 1007;
    else if (std::fabs(alpha - 11.97343036) < epsilon)  return 2031;
    else if (std::fabs(alpha - 12.9613967) < epsilon)   return 4031;
    else if (std::fabs(alpha - 13.96358391) < epsilon)  return 8127;
    else return 0;  // Default case when alpha does not match
}


usint alpha2depthOpenFHE(const float alpha)
{
    const float epsilon = 0.0001f; // Tolerance for comparing floating-point values
    
    if (std::fabs(alpha - 3.091503713) < epsilon)       return 3;
    else if (std::fabs(alpha - 3.62961696) < epsilon)   return 3;
    else if (std::fabs(alpha - 4.026158059) < epsilon)  return 4;
    else if (std::fabs(alpha - 4.339009337) < epsilon)  return 4;
    else if (std::fabs(alpha - 4.596911194) < epsilon)  return 4;
    else if (std::fabs(alpha - 4.816362316) < epsilon)  return 4;
    else if (std::fabs(alpha - 5.006872028) < epsilon)  return 5;
    else if (std::fabs(alpha - 5.175246808) < epsilon)  return 5;
    else if (std::fabs(alpha - 5.326378453) < epsilon)  return 5;
    else if (std::fabs(alpha - 5.449221655) < epsilon)  return 5;
    else if (std::fabs(alpha - 5.575712623) < epsilon)  return 5;
    else if (std::fabs(alpha - 5.69235669) < epsilon)   return 5;
    else if (std::fabs(alpha - 5.800459805) < epsilon)  return 5;
    else if (std::fabs(alpha - 5.900865483) < epsilon)  return 6;
    else if (std::fabs(alpha - 5.99483256) < epsilon)   return 6;
    else if (std::fabs(alpha - 6.082714136) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.166037129) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.244093359) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.318384415) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.389465946) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.448095937) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.510721532) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.573560838) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.63309892) < epsilon)   return 6;
    else if (std::fabs(alpha - 6.690314796) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.745604685) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.799211286) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.849537749) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.899942561) < epsilon)  return 6;
    else if (std::fabs(alpha - 6.946884719) < epsilon)  return 7;
    else if (std::fabs(alpha - 6.993709998) < epsilon)  return 7;
    else if (std::fabs(alpha - 8.939735521) < epsilon)  return 8;
    else if (std::fabs(alpha - 9.939149282) < epsilon)  return 9;
    else if (std::fabs(alpha - 10.962565805) < epsilon)  return 10;
    else if (std::fabs(alpha - 11.97343036) < epsilon)  return 11;
    else if (std::fabs(alpha - 12.9613967) < epsilon)   return 12;
    else if (std::fabs(alpha - 13.96358391) < epsilon)  return 14;
    else return 0;  // Default case when alpha does not match
}

usint degree2depthOpenFHE(const int degree) {
    switch(degree) {
        case 3: return 3;
        case 5: return 3;
        case 11: return 4;
        case 12: return 4;
        case 13: return 4;
        case 15: return 5;
        case 17: return 5;
        case 21: return 5;
        case 25: return 5;
        case 51: return 6;
        case 101: return 7;
        case 111: return 7;
        case 125: return 8;
        case 143: return 8;
        case 201: return 8;
        case 251: return 9;
        case 401: return 9;
        case 501: return 10;
        case 1001: return 10;
        case 1113: return 11;
        case 1601: return 11;
        case 1671: return 11;
        case 2007: return 11;
        case 3201: return 12;
        case 3355: return 12;
        case 6401: return 13;
        case 12801: return 14;
        default: return 0;
    }
}

float degree2errorOpenFHE(const int degree){
    switch(degree) {
        case 3: return 0.3;
        case 5: return 0.2;
        case 11: return 0.1;
        case 12: return 0.09;
        case 13: return 0.08;
        case 15: return 0.07;
        case 17: return 0.06;
        case 21: return 0.05;
        case 25: return 0.06;
        case 51: return 0.02;
        case 101: return 0.01;
        case 111: return 0.009;
        case 125: return 0.008;
        case 143: return 0.007;
        case 201: return 0.005;
        case 251: return 0.004;
        case 401: return 0.003;
        case 501: return 0.002;
        case 1001: return 0.001;
        case 1113: return 0.0009;
	case 1601: return 0.0007;
        case 1671: return 0.0006;
        case 2007: return 0.0005;
        case 3201: return 0.0004;
        case 3355: return 0.0003;
        case 6401: return 0.0002;
        case 12801: return 0.0001;
        default: return 0;
    }
}
std::vector<int> alpha2degreeCompos(const float alpha) {

    // Use epsilon to compare floating-point numbers
    const float epsilon = 0.0001;


    if (std::fabs(alpha - 2.000) < epsilon) { 
        return {3};
    } else if (std::fabs(alpha - 2.321928095) < epsilon) {
        return {3};
    } else if (std::fabs(alpha - 2.736965594) < epsilon) {
        return {7};
    } else if (std::fabs(alpha - 3.321928095) < epsilon) {
        return {15};
    } else if (std::fabs(alpha - 4.321928095) < epsilon) {
        return {7,7};
    } else if (std::fabs(alpha - 4.473931188) < epsilon) {
        return {7,7};
    } else if (std::fabs(alpha - 4.64385619) < epsilon) {
        return {47};
        //return {7,7};
    } else if (std::fabs(alpha - 4.836501268) < epsilon) {
        return {7,11};
    } else if (std::fabs(alpha - 5.058893689) < epsilon) {
        return {7,13};
    } else if (std::fabs(alpha - 5.321928095) < epsilon) {
        return {3, 7, 7};
    } else if (std::fabs(alpha - 5.64385619) < epsilon) {
        return {13, 13};
    } else if (std::fabs(alpha - 6.058893689) < epsilon) {
        return {15, 15};
    } else if (std::fabs(alpha - 6.64385619) < epsilon) {
        return {5, 7, 13};
    } else if (std::fabs(alpha - 7.64385619) < epsilon) {
        return {3, 5, 7, 13};
    } else if (std::fabs(alpha - 7.795859283) < epsilon) {
        return {13, 7, 15};
    } else if (std::fabs(alpha - 7.965784285) < epsilon) {
        return {15, 7, 15};
    } else if (std::fabs(alpha - 8.158429363) < epsilon) {
        return {15, 7, 3, 7};
    } else if (std::fabs(alpha - 8.380821784) < epsilon) {
        return {15, 11, 13};
    } else if (std::fabs(alpha - 8.64385619) < epsilon) {
        return {15, 13, 15};
    } else if (std::fabs(alpha - 8.965784285) < epsilon) {
        return {15, 3, 7, 15};
    } else if (std::fabs(alpha - 9.380821784) < epsilon) {
        return {15, 15, 27};
    } else if (std::fabs(alpha - 9.965784285) < epsilon) {
        return {15, 19, 27};
    } else if (std::fabs(alpha - 10.96578428) < epsilon) {
        return {15, 7, 9, 23};
    } else if (std::fabs(alpha - 11.11778738) < epsilon) {
        return {15, 9, 7, 23};
    } else if (std::fabs(alpha - 11.28771238) < epsilon) {
        return {15, 9, 9, 23};
    } else if (std::fabs(alpha - 11.48035746) < epsilon) {
        return {15, 9, 9, 27};
    } else if (std::fabs(alpha - 11.70274988) < epsilon) {
        return {15, 9, 9, 31};
    } else if (std::fabs(alpha - 11.96578428) < epsilon) {
        return {15, 27, 55};
    } else if (std::fabs(alpha - 12.28771238) < epsilon) {
        return {15, 31, 61};
    } else if (std::fabs(alpha - 12.70274988) < epsilon) {
        return {15, 59, 59};
    } else if (std::fabs(alpha - 13.28771238) < epsilon) {
        return {15, 15, 11, 27};
    } else if (std::fabs(alpha - 14.28771238) < epsilon) {
        return {15, 15, 27, 27};
    } else {
        return {};  // Return an empty vector if alpha does not match
    }
}


std::vector<int> alpha2degreeComposOld(const float alpha) {

    // Use epsilon to compare floating-point numbers
    const float epsilon = 0.0001;

    if (std::fabs(alpha - 2.000) < epsilon) { return {3};
    } else if (std::fabs(alpha - 2.321928095) < epsilon) {
        return {3};
    } else if (std::fabs(alpha - 2.736965594) < epsilon) {
        return {3, 3};
    } else if (std::fabs(alpha - 3.321928095) < epsilon) {
        return {3, 5};
    } else if (std::fabs(alpha - 4.321928095) < epsilon) {
        return {3, 3, 3, 3};
    } else if (std::fabs(alpha - 4.473931188) < epsilon) {
        return {3, 3, 3, 3};
    } else if (std::fabs(alpha - 4.64385619) < epsilon) {
        return {3, 3, 3, 3};
        //return {7,7};
    } else if (std::fabs(alpha - 4.836501268) < epsilon) {
        return {3, 3, 3, 5};
    } else if (std::fabs(alpha - 5.058893689) < epsilon) {
        return {3, 3, 3, 5};
    } else if (std::fabs(alpha - 5.321928095) < epsilon) {
        return {3, 3, 3, 3, 3};
    } else if (std::fabs(alpha - 5.64385619) < epsilon) {
        return {3, 3, 3, 3, 3};
    } else if (std::fabs(alpha - 6.058893689) < epsilon) {
        return {3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 6.64385619) < epsilon) {
        return {3, 3, 3, 3, 3, 3};
    } else if (std::fabs(alpha - 7.64385619) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3};
    } else if (std::fabs(alpha - 7.795859283) < epsilon) {
        return {3, 3, 3, 3, 5, 5};
    } else if (std::fabs(alpha - 7.965784285) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 8.158429363) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 8.380821784) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 3};
    } else if (std::fabs(alpha - 8.64385619) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 3};
    } else if (std::fabs(alpha - 8.965784285) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 9.380821784) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 9.965784285) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 10.96578428) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    } else if (std::fabs(alpha - 11.11778738) < epsilon) {
       return {3, 3, 3, 3, 3, 3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 11.28771238) < epsilon) {
       return {3, 3, 3, 3, 3, 3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 11.48035746) < epsilon) {
       return {3, 3, 3, 3, 3, 3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 11.70274988) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 11.96578428) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    } else if (std::fabs(alpha - 12.28771238) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 12.70274988) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 3, 3, 5};
    } else if (std::fabs(alpha - 13.28771238) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    } else if (std::fabs(alpha - 14.28771238) < epsilon) {
        return {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    } else {
        std::cout << "Alpha value not found." << std::endl;
        return {};
    }
}


float alpha2error(const float alpha) {
    
    // Use epsilon to compare floating-point numbers
    const float epsilon = 0.0001;

    if (std::fabs(alpha - 2.000) < epsilon) { return {0.5};
    } else if (std::fabs(alpha - 2.321928095) < epsilon) {
        return 0.4;
    } else if (std::fabs(alpha - 2.736965594) < epsilon) {
        return 0.3;
    } else if (std::fabs(alpha - 3.321928095) < epsilon) {
        return 0.2;
    } else if (std::fabs(alpha - 4.321928095) < epsilon) {
        return 0.1;
    } else if (std::fabs(alpha - 4.473931188) < epsilon) {
        return 0.09;
    } else if (std::fabs(alpha - 4.64385619) < epsilon) {
        return 0.08;
    } else if (std::fabs(alpha - 4.836501268) < epsilon) {
        return 0.07;
    } else if (std::fabs(alpha - 5.058893689) < epsilon) {
        return 0.06;
    } else if (std::fabs(alpha - 5.321928095) < epsilon) {
        return 0.05;
    } else if (std::fabs(alpha - 5.64385619) < epsilon) {
        return 0.04;
    } else if (std::fabs(alpha - 6.058893689) < epsilon) {
        return 0.03;
    } else if (std::fabs(alpha - 6.64385619) < epsilon) {
        return 0.02;
    } else if (std::fabs(alpha - 7.64385619) < epsilon) {
        return 0.01;
    } else if (std::fabs(alpha - 7.795859283) < epsilon) {
        return 0.009;
    } else if (std::fabs(alpha - 7.965784285) < epsilon) {
        return 0.008;
    } else if (std::fabs(alpha - 8.158429363) < epsilon) {
        return 0.007;
    } else if (std::fabs(alpha - 8.380821784) < epsilon) {
        return 0.006;
    } else if (std::fabs(alpha - 8.64385619) < epsilon) {
        return 0.005;
    } else if (std::fabs(alpha - 8.965784285) < epsilon) {
        return 0.004;
    } else if (std::fabs(alpha - 9.380821784) < epsilon) {
        return 0.003;
    } else if (std::fabs(alpha - 9.965784285) < epsilon) {
        return 0.002;
    } else if (std::fabs(alpha - 10.96578428) < epsilon) {
        return 0.001;
    } else if (std::fabs(alpha - 11.11778738) < epsilon) {
        return 0.0009;
    } else if (std::fabs(alpha - 11.28771238) < epsilon) {
        return 0.0008;
    } else if (std::fabs(alpha - 11.48035746) < epsilon) {
        return 0.0007;
    } else if (std::fabs(alpha - 11.70274988) < epsilon) {
        return 0.0006;
    } else if (std::fabs(alpha - 11.96578428) < epsilon) {
        return 0.0005;
    } else if (std::fabs(alpha - 12.28771238) < epsilon) {
        return 0.0004;
    } else if (std::fabs(alpha - 12.70274988) < epsilon) {
        return 0.0003;
    } else if (std::fabs(alpha - 13.28771238) < epsilon) {
        return 0.0002;
    } else if (std::fabs(alpha - 14.28771238) < epsilon) {
        return 0.0001;
    } else {
        std::cout << "Alpha value not found." << std::endl;
        return {};
    }

}


usint alpha2depthCompos(const float alpha)
{
    // Use epsilon to compare floating-point numbers
    const float epsilon = 0.0001;

    if (std::fabs(alpha - 2.000) < epsilon) { 
        return 5;
    } else if (std::fabs(alpha - 2.321928095) < epsilon) {
        return 5;
    } else if (std::fabs(alpha - 2.736965594) < epsilon) {
        return 5;
    } else if (std::fabs(alpha - 3.321928095) < epsilon) {
        return 5;
    } else if (std::fabs(alpha - 4.321928095) < epsilon) {
        return 8;
    } else if (std::fabs(alpha - 4.473931188) < epsilon) {
        return 8;
    } else if (std::fabs(alpha - 4.64385619) < epsilon) {
        return 5;
    } else if (std::fabs(alpha - 4.836501268) < epsilon) {
        return 8;
    } else if (std::fabs(alpha - 5.058893689) < epsilon) {
        return 8;
    } else if (std::fabs(alpha - 5.321928095) < epsilon) {
        return 11;
    } else if (std::fabs(alpha - 5.64385619) < epsilon) {
        return 10;
    } else if (std::fabs(alpha - 6.058893689) < epsilon) {
        return 8;
    } else if (std::fabs(alpha - 6.64385619) < epsilon) {
        return 12;
    } else if (std::fabs(alpha - 7.64385619) < epsilon) {
        return 15;
    } else if (std::fabs(alpha - 7.795859283) < epsilon) {
        return 12;
    } else if (std::fabs(alpha - 7.965784285) < epsilon) {
        return 12;
    } else if (std::fabs(alpha - 8.158429363) < epsilon) {
        return 16;
    } else if (std::fabs(alpha - 8.380821784) < epsilon) {
        return 13;
    } else if (std::fabs(alpha - 8.64385619) < epsilon) {
        return 13;
    } else if (std::fabs(alpha - 8.965784285) < epsilon) {
        return 16;
    } else if (std::fabs(alpha - 9.380821784) < epsilon) {
        return 12;
    } else if (std::fabs(alpha - 9.965784285) < epsilon) {
        return 12;
    } else if (std::fabs(alpha - 10.96578428) < epsilon) {
        return 16;
    } else if (std::fabs(alpha - 11.11778738) < epsilon) {
        return 16;
    } else if (std::fabs(alpha - 11.28771238) < epsilon) {
        return 16;
    } else if (std::fabs(alpha - 11.48035746) < epsilon) {
        return 17;
    } else if (std::fabs(alpha - 11.70274988) < epsilon) {
        return 17;
    } else if (std::fabs(alpha - 11.96578428) < epsilon) {
        return 14;
    } else if (std::fabs(alpha - 12.28771238) < epsilon) {
        return 14;
    } else if (std::fabs(alpha - 12.70274988) < epsilon) {
        return 19;
    } else if (std::fabs(alpha - 13.28771238) < epsilon) {
        return 17;
    } else if (std::fabs(alpha - 14.28771238) < epsilon) {
        return 17;
    } else {
        std::cout << "Alpha value not found." << std::endl;
        return 5;
    }

}


std::vector<double> read_coefficients_string(
    const std::string &coefficient_path
)
{
    // Read float values from the file
    std::ifstream infile(coefficient_path);
    std::vector<float> coefficients;
    float value;
    while (infile >> value) {
        coefficients.push_back(value);
    }
    infile.close();

    // Convert float values to double
    //std::vector<double> double_coefficients(coefficients.begin(), coefficients.end());
    std::vector<double> double_coefficients;
    double_coefficients.reserve(coefficients.size()); // Reserve space for efficiency
    for (float value : coefficients) {
            double_coefficients.push_back(static_cast<double>(value));
    }

    return double_coefficients;
}


Ciphertext<DCRTPoly> min_func(
    const Ciphertext<DCRTPoly> &c1,
    const Ciphertext<DCRTPoly> &c2,
    double a,
    double b,
    uint32_t degree,
    std::vector<std::vector<double>> coefficient_list,
    std::vector<int> degrees,
    double error
)
{
    //auto add = c1 + c2;
    auto add = c1->GetCryptoContext()->EvalAdd(c1, c2);
    //auto sub = c1 - c2;
    auto sub = c1->GetCryptoContext()->EvalSub(c1, c2);

    auto sign = sub;
    for(size_t i = 0; i < degrees.size(); i++) {
        std::vector<double> double_coefficients = coefficient_list[i];

        sign = c1->GetCryptoContext()->EvalChebyshevSeries(sign, double_coefficients, a, b);
    }
    
    auto mult = c1->GetCryptoContext()->EvalMult(sign, sub);
    auto max = c1->GetCryptoContext()->EvalAdd(add, mult);
    max = c1->GetCryptoContext()->EvalMult(max, 0.5);
    auto min = c1->GetCryptoContext()->EvalSub(add, max);
    return min;
    // return sign;
}



Ciphertext<DCRTPoly> compare_composite(
    const Ciphertext<DCRTPoly> &c1,
    const Ciphertext<DCRTPoly> &c2,
    std::vector<std::vector<double>> coefficient_list,
    double a,
    double b,
    std::vector<int> M_degs,
    float alpha,
    std::string alpha_string,
    double error
)
{

    auto resultC = c1 - c2;
    std::string output_directory = "/home/s2086034/remez_outputs";;
    for (size_t q = 0; q < M_degs.size(); ++q) {
        std::vector<double> coefficients = coefficient_list[q];
        resultC = c1->GetCryptoContext()->EvalChebyshevSeries(resultC, coefficients, -1, 1);
    }

    resultC = resultC + 1;

    resultC = c1->GetCryptoContext()->EvalMult(resultC, 0.5);

    return resultC;

}



Ciphertext<DCRTPoly> compare(
    const Ciphertext<DCRTPoly> &c1,
    const Ciphertext<DCRTPoly> &c2,
    double a,
    double b,
    uint32_t degree,
    double error
)
{
    return c1->GetCryptoContext()->EvalChebyshevFunction(
        [error](double x) -> double { 
            if      (x > error)   return 1;
            else if (x >= -error) return 0.5;
            else                  return 0;
        },
        c1 - c2,
        a, b, degree
    );
}



Ciphertext<DCRTPoly> compareGt(
    const Ciphertext<DCRTPoly> &c1,
    const Ciphertext<DCRTPoly> &c2,
    double a,
    double b,
    uint32_t degree,
    double error
)
{
    return c1->GetCryptoContext()->EvalChebyshevFunction(
        [error](double x) -> double { 
            if      (x > error)   return 1;
            else                  return 0;
        },
        c1 - c2,
        a, b, degree
    );
}


Ciphertext<DCRTPoly> indicator(
    const Ciphertext<DCRTPoly> &c,
    double a1,
    double b1,
    double a,
    double b,
    uint32_t degree
)
{
    return c->GetCryptoContext()->EvalChebyshevFunction(
        [a1,b1](double x) -> double {
            return (x < a1 || x > b1) ? 0 : 1; },
        c,
        a, b, degree
    );
}


// int main()
// {

//     const usint vectorLength            = 4;
//     const usint functionDepth           = 10;

//     const usint integralPrecision       = 10;
//     const usint decimalPrecision        = 50;
//     const usint multiplicativeDepth     = functionDepth;
//     const usint numSlots                = vectorLength;
//     const bool enableBootstrap          = false;
//     const usint ringDim                 = 0;
//     const bool verbose                  = true;

//     std::vector<int32_t> indices = {};

//     CryptoContext<DCRTPoly> cryptoContext = generateCryptoContext(
//         integralPrecision,
//         decimalPrecision,
//         multiplicativeDepth,
//         numSlots,
//         enableBootstrap,
//         ringDim,
//         verbose
//     );

//     KeyPair<DCRTPoly> keyPair = keyGeneration(
//         cryptoContext,
//         indices,
//         numSlots,
//         enableBootstrap,
//         verbose
//     );

//     std::vector<double> v = {-0.75, -0.1, 0.0, 0.75};
//     std::vector<double> zero = {0.0, 0.0, 0.0, 0.0};

//     std::cout << "Vector: " << v << std::endl;

//     Ciphertext<DCRTPoly> vC = cryptoContext->Encrypt(
//         keyPair.publicKey,
//         cryptoContext->MakeCKKSPackedPlaintext(v)
//     );
//     Ciphertext<DCRTPoly> zeroC = cryptoContext->Encrypt(
//         keyPair.publicKey,
//         cryptoContext->MakeCKKSPackedPlaintext(zero)
//     );

//     auto start = std::chrono::high_resolution_clock::now();

//     // Ciphertext<DCRTPoly> resultC = compare(
//     //     zeroC, vC,
//     //     -1.0, 1.0,
//     //     depth2degree(functionDepth)
//     // );
//     Ciphertext<DCRTPoly> resultC = indicator(
//         vC,
//         -0.5, 0.5,
//         -1.0, 1.0,
//         depth2degree(functionDepth)
//     );

//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_seconds = end - start;
//     std::cout << elapsed_seconds.count() << "s" << std::endl;

//     Plaintext resultP;
//     cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
//     resultP->SetLength(vectorLength);

//     std::vector<double> result = resultP->GetRealPackedValue();
//     std::cout << "Result: " << result << std::endl;

//     return 0;

// }
