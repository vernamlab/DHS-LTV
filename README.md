# LaSH (Leveled and Somehat Homomorphic) Encryption Library

The library provides a homomorphic encryption scheme by Doroz, Hu and Sunar (DHS), along with several applications (coming soon). The details of the scheme can be seen in : https://eprint.iacr.org/2014/039.pdf. Our scheme is a single key variant of Lopez-Alt, Tromer and Vaikuntanathan (LTV) Homomorphic Encryption. The original construction can be found in : https://eprint.iacr.org/2013/094.pdf.

# Somewhat and leveled Fully Homomorphic Encryption

Our library implements SHE and leveled-FHE versions of the same scheme. We will provide a setup tutorial soon. This choice is mostly important for application design. Depending on the depth of the homomorphic circuit that we want to evaluate, we can make a decision between the two, for a more efficient design.

# Applications

Using the core homomorphic cryptosystem, we have built several applications. Examples include encrypted computation of cryptographic protocols such as AES, Prince Cipher, as well as more practical applications like PIR, Sorting Encrypted Data, Arithmetic over Encrypted Data and encrypted keyword search/completion. Sorting application is already available can be seen in applications folder, the other source codes will be avaliabe soon. All of these projects were developed for research purpose only, therefore although they provide maximum efficiency, some parts may lack a clean user interface for now. A tutorial document is in progress. The links for the mentioned applications:
  1. PIR      :https://eprint.iacr.org/2014/232.pdf
  2. PRINCE   :https://eprint.iacr.org/2014/233.pdf
  3. SORTING  :https://eprint.iacr.org/2015/274.pdf
  4. AUTOCOMPLETE :https://eprint.iacr.org/2015/1194.pdf

# Installation

This library uses NTL (v9.0.2 or later) with GMP support (with C++11 on). Tested on Linux environment. There is a make file available.



