# Protparam
# ProtParam - C++ Implementation

A lightweight C++ implementation of the **ExPASy ProtParam** web tool for analyzing protein sequences. This project calculates several important physicochemical properties of proteins using only their amino acid sequence.

## Overview

ProtParam is a popular bioinformatics tool used to estimate various biochemical properties of proteins. This project recreates some of its core functionalities in C++, allowing offline analysis without requiring access to the ExPASy web server.

## Features

- Molecular Weight calculation
- Amino Acid Composition
- Count of Acidic and Basic Amino Acids
- Theoretical Isoelectric Point (pI) *(simplified implementation)*
- Extinction Coefficient estimation
- Fast execution with minimal dependencies
- Console-based interface

## Implemented Properties

| Property | Status |
|----------|--------|
| Molecular Weight | ✅ |
| Amino Acid Composition | ✅ |
| Acidic/Basic Residue Count | ✅ |
| Theoretical pI | ✅ (Simplified) |
| Extinction Coefficient | ✅ |
| Instability Index | ❌ |
| Aliphatic Index | ❌ |
| GRAVY | ❌ |
| Estimated Half-life | ❌ |

## Project Structure

```
Protparam-project/
│
├── main.cpp          # Main program
├── Protparam.cpp     # Property calculation functions
└── README.md
```

## Requirements

- C++11 or newer
- GCC / Clang / MSVC

## Compilation

Using g++:

```bash
g++ main.cpp Protparam.cpp -o protparam
```

or

```bash
g++ *.cpp -o protparam
```

## Usage

Run the executable:

```bash
./protparam
```

Example input:

```
Enter the protein sequence:
MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE
```

Example output:

```
Molecular Weight: 5123.54 Da

Amino Acid Composition
A: 2
D: 1
F: 5
...

Number of acidic amino acids: 3
Number of basic amino acids: 8

Theoretical pI: 8.67

Extinction Coefficient: 12480
```

## Algorithms Used

- Amino acid molecular weight lookup table
- Frequency counting using `std::unordered_map`
- Simplified pKa-based pI estimation
- Extinction coefficient calculation based on:
  - Tryptophan (W)
  - Tyrosine (Y)
  - Cysteine (C)

## Limitations

This project is intended for educational and learning purposes.

Current limitations include:

- Simplified pI calculation
- No support for FASTA files
- No batch sequence processing
- Missing Instability Index
- Missing Aliphatic Index
- Missing GRAVY calculation
- No graphical interface

## Future Work

- Add FASTA file support
- Implement the complete ProtParam algorithm
- Add Instability Index
- Add Aliphatic Index
- Add GRAVY calculation
- Improve theoretical pI accuracy
- Export results to CSV
- Unit testing
- CMake support

## Disclaimer

This project is an independent implementation inspired by the **ExPASy ProtParam** web application. It is intended for educational and research purposes only and is not affiliated with or endorsed by ExPASy.

## Author

**Praveen Gautam**

GitHub: https://github.com/Praveengautam759

## License

This project is released under the MIT License.
