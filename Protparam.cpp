#include <iostream>
#include <string>
#include <unordered_map>
#include <iomanip>
#include <cmath>

using namespace std;
// Molecular weights of amino acids
unordered_map<char, double> amino_acid_weights = {
    {'A', 89.09}, {'R', 174.20}, {'N', 132.12}, {'D', 133.10},
    {'C', 121.15}, {'E', 147.13}, {'Q', 146.15}, {'G', 75.07},
    {'H', 155.16}, {'I', 131.17}, {'L', 131.17}, {'K', 146.19},
    {'M', 149.21}, {'F', 165.19}, {'P', 115.13}, {'S', 105.09},
    {'T', 119.12}, {'W', 204.23}, {'Y', 181.19}, {'V', 117.15}
};

// pKa values for acidic and basic amino acids (simplified)
unordered_map<char, double> pKa_values = {
    {'D', 3.9},  {'E', 4.2},  // Acidic amino acids
    {'C', 8.3},  {'H', 6.0},  // Weakly acidic
    {'K', 10.5}, {'R', 12.5}  // Basic amino acids
};

// Extinction coefficients for W, Y, and disulfide bonds (C)
unordered_map<char, double> extinction_coefficients = {
    {'W', 5500}, {'Y', 1490}, {'C', 125}
};

// Function to calculate molecular weight
double calculateMolecularWeight(const string& sequence) {
    double totalWeight = 0.0;
    for (char aa : sequence) {
        totalWeight += amino_acid_weights[aa];
    }
    return totalWeight;
}

// Function to calculate amino acid composition
unordered_map<char, int> calculateAAComposition(const string& sequence) {
    std::unordered_map<char, int> composition;
    for (char aa : sequence) {
        composition[aa]++;
    }
    return composition;
}

// Function to calculate the number of acidic and basic amino acids
pair<int, int> countAcidicBasic(const unordered_map<char, int>& composition) {
    int acidic = composition.at('D') + composition.at('E');
    int basic = composition.at('K') + composition.at('R') + composition.at('H');
    return {acidic, basic};
}

// Function to calculate extinction coefficient
double calculateExtinctionCoefficient(const unordered_map<char, int>& composition) {
    return (composition.at('W') * extinction_coefficients['W']) +
           (composition.at('Y') * extinction_coefficients['Y']) +
           (composition.at('C') * extinction_coefficients['C']);
}

// Function to calculate theoretical pI (simplified)
double calculatePI(const unordered_map<char, int>& composition) {
    double pH = 7.0;
    double charge = 0.0;
    while (pH <= 14.0) {
        charge = 0.0;
        for (auto& pKa : pKa_values) {
            if (composition.find(pKa.first) != composition.end()) {
                if (pKa.first == 'D' || pKa.first == 'E') {
                    charge -= composition.at(pKa.first) / (1 + pow(10, pH - pKa.second)); // acidic groups
                } else {
                    charge += composition.at(pKa.first) / (1 + pow(10, pKa.second - pH)); // basic groups
                }
            }
        }
        if (fabs(charge) < 0.01) break;
        pH += 0.01;
    }
    return pH;
}

int main() {
    string protein_sequence;
    cout << "Enter the protein sequence (use single letter amino acid codes): ";
    cin >> protein_sequence;

    // Molecular weight
    double molecular_weight = calculateMolecularWeight(protein_sequence);
    cout << "Molecular Weight: " << fixed << setprecision(2) << molecular_weight << " Da" << endl;

    // Amino acid composition
    auto composition = calculateAAComposition(protein_sequence);
    cout << "Amino Acid Composition:" << endl;
    for (auto& aa : composition) {
        cout << aa.first << ": " << aa.second << endl;
    }

    // Number of acidic and basic amino acids
    auto [acidic, basic] = countAcidicBasic(composition);
    cout << "Number of acidic amino acids: " << acidic << endl;
    cout << "Number of basic amino acids: " << basic << endl;

    // Theoretical pI
    double theoretical_pI = calculatePI(composition);
    cout << "Theoretical pI: " << fixed << setprecision(2) << theoretical_pI << endl;

    // Extinction coefficient
    double extinction_coeff = calculateExtinctionCoefficient(composition);
    cout << "Extinction Coefficient: " << extinction_coeff << endl;

    return 0;
}
