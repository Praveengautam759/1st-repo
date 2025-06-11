#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <cctype>

using namespace std;

// Amino acid molecular weights (average, in Daltons)
map<char, double> aminoAcidMass = {
    {'A',  89.09}, {'R', 174.20}, {'N', 132.12}, {'D', 133.10},
    {'C', 121.15}, {'Q', 146.15}, {'E', 147.13}, {'G',  75.07},
    {'H', 155.16}, {'I', 131.17}, {'L', 131.17}, {'K', 146.19},
    {'M', 149.21}, {'F', 165.19}, {'P', 115.13}, {'S', 105.09},
    {'T', 119.12}, {'W', 204.23}, {'Y', 181.19}, {'V', 117.15}
};

// Extinction coefficients for W, Y, C (in M^-1 cm^-1)
const int extW = 5500;
const int extY = 1490;
const int extC = 125; // for disulfide (cystine)

bool isValidAminoAcid(char aa) {
    return aminoAcidMass.find(aa) != aminoAcidMass.end();
}

int main() {
    string rawInput, sequence;
    map<char, int> aaCount;

    cout << "Enter the protein sequence (single-letter amino acid codes):\n";
    getline(cin, rawInput);

    // Clean and validate sequence
    for (char ch : rawInput) {
        ch = toupper(ch);
        if (isValidAminoAcid(ch)) {
            sequence += ch;
            aaCount[ch]++;
        }
    }

    if (sequence.empty()) {
        cerr << "Invalid or empty sequence. Only A-Z amino acids are allowed." << endl;
        return 1;
    }

    int length = sequence.length();
    double molWeight = 0.0;
    int acidicCount = aaCount['D'] + aaCount['E'];
    int basicCount = aaCount['K'] + aaCount['R'];
    int cCount = aaCount['C'], yCount = aaCount['Y'], wCount = aaCount['W'];

    // Molecular weight calculation
    for (auto& pair : aaCount) {
        molWeight += aminoAcidMass[pair.first] * pair.second;
    }
    // Subtract water (18.015 Da) per peptide bond (n-1)
    molWeight -= (length - 1) * 18.015;

    // Extinction coefficients
    int extCoeffReduced = wCount * extW + yCount * extY;
    int extCoeffCystine = extCoeffReduced + (cCount / 2) * extC; // assuming all Cys form pairs

    // Output
    cout << "\n========== Protein Analysis ==========" << endl;
    cout << "Protein Length: " << length << " residues" << endl;
    cout << fixed << setprecision(2);
    cout << "Molecular Weight: " << molWeight << " Da" << endl;
    cout << "Extinction Coefficient (reduced Cys): " << extCoeffReduced << " M⁻¹cm⁻¹" << endl;
    cout << "Extinction Coefficient (disulfide/Cystine): " << extCoeffCystine << " M⁻¹cm⁻¹" << endl;

    double acidicPercent = (double)acidicCount / length * 100;
    double basicPercent = (double)basicCount / length * 100;
    cout << "Percentage of acidic residues (D+E): " << acidicPercent << "%" << endl;
    cout << "Percentage of basic residues (K+R): " << basicPercent << "%" << endl;

    cout << "\nAmino Acid Composition:\n";
    for (auto& pair : aaCount) {
        double percent = (double)pair.second / length * 100;
        cout << "  " << pair.first << ": " << pair.second << " (" << percent << "%)" << endl;
    }

    return 0;
}
