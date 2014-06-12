// rsa.cpp
// CMPSC101 Final Project
// Ken Hwang (krh5058)

#include <algorithm>
#include <locale>
#include <iostream>
using namespace std;

#pragma warning(disable: 4800)
#include <mpirxx.h>
#pragma warning(default: 4800)

#include <sstream>
#include <string>
#include <time.h> 
#include <vector>

void primes(vector<int> &in_vector, int n){
	/*
		Get all primes up to integer n.
		Check against previous primes, and add to input vector if more primes are found
		Assume input vector includes first prime number 2.
	*/
    for(int i=3; i < n; i+=2) // Iterate by odds
    {
        bool isPrime=true;
        for(int p=0; p<in_vector.size(); p++) // Index through previous primes
        {
            if(i % in_vector[p] == 0) // If divisible, then not prime
            {
                isPrime=false;
                break;
            }
        }
        if(isPrime) 
        {
            in_vector.push_back(i);
        }
    }
}

vector<int> factors(vector<int> &n_primes, int n){
	/*
		Get all factors for integer n.
		Iterate through vector of primes for integer n, and check if division results in an integer
		Output factor_vector
	*/
	vector<int> factor_vector;
	for(int i=0; i<n_primes.size(); i++)
	{
		if(n % n_primes[i] == 0)
		{
			factor_vector.push_back(n_primes[i]);
		}
	}
	return factor_vector;
}

void get_alphabet(vector<char> &alphabet){
	/*
		Fill vector<char> with uppercase letters in alphabet
		Iterate through ASCII codes 65 to 90
	*/
	int a=65; // Start with ASCII for "A"
	while(a<=90) // End with ASCII for "Z"
	{
		alphabet.push_back((char)a);
		//cout << a << " " << (char)a << endl;
		a++;
	}
}

bool check_rel_prime(int a, int b){
	/*
		Check if integers "a" and "b" are relatively prime
		Iterate, starting at 2, until "i" becomes larger than one of the values
		Check against either value to determine if they are divisible by "i"
		If divisible, then these numbers are not relatively prime
	*/
	bool rel_prime = true;
	for(int i=2; i <= a && i <= b; i++)
	{
		if(a % i == 0 && b % i == 0)
		{
			//cout << "not rel prime" << endl;
			rel_prime = false;
			break;
		}
	}
	return rel_prime;
}

int main (void)
{
	int val = 0;

	// Intro
	cout << "RSA is an encryption algorithm that uses a public key and a private key to encode and decode messages, respectively." << endl;
	cout << "We will demonstrate how this works in a very basic fashion." << endl;
	cout << "To start, we must first generate these keys." << endl;
	cout << "Press Enter to Continue..." << endl;
	cin.ignore();
	
	// Query user
	cout << "We need to determine two starting values, U and V." << endl;
	cout << "To make it difficult to decypher, we want these values to be prime numbers." << endl;
	cout << "To begin, let's specify an upper range of prime numbers and randomly pick two values from that range." << endl;
	cout << "Enter an integer between 20 and 300." << endl;
	while (val < 20 || val > 300)
	{
		cin >> val;
		if(val < 20 || val > 300)
		{
			cout << "Enter an integer greater than 19 and less than 301." << endl;
		}
	}

	cout << endl;
	cin.ignore();

	// Generate and display prime number vector
	cout << "Great! You chose " << val << "." << endl;
	cout << "Generating a prime number vector between 2 and " << val << "..." << endl;
    vector<int> prime_vector;
    prime_vector.push_back(2); // Enter 2
    prime_vector.push_back(3); // Enter 3
	primes(prime_vector,val);
	
	cout << "The prime numbers we can choose from are..." << endl;
	for(int i=0; i < prime_vector.size(); i++)
	{
		cout << prime_vector[i] << " ";
	}
	cout << endl;
	cout << endl;

	// Pick two prime numbers and display
	cout << "Let's randomly pick 2 prime numbers." << endl;	
	cout << "Press Enter to Continue..." << endl;
	cin.ignore();

	srand (time(NULL)); // Random Seed
	
	int index1 = rand() % (prime_vector.size()-3) + 3; // Primes must be 7 (prime_vector[3]) or larger,
	int index2 = rand() % (prime_vector.size()-3) + 3; // Primes must be 7 (prime_vector[3]) or larger
	while (index1 == index2){ // If index2 is the same value
		index2 = rand() % (prime_vector.size()-3) + 3;
	}
	int U = prime_vector[index1];
	int V = prime_vector[index2];
	cout << "We have chosen " << U << " and " << V << "." << endl;
	cout << "These values will be U and V, respectively." << endl;	
	cout << "Next, we can use these two prime numbers to determine R, the modulo." << endl;
	cout << "All this takes is to multiply U and V (U*V=R)." << endl;
	cout << "Press Enter to Continue..." << endl;
	cin.ignore();

	// Calculate R and phi(R)
	int R = U*V;

	cout << U << " times " << V << " equals " << R << endl;
	cout << endl;
	cout << "This seems simple, but the mathematics behind RSA don't only involve R, but phi(R)." << endl;
	cout << "phi() is Euler's totient function, and is important to how RSA works." << endl;
	cout << "In brief, phi(N) is the total amount of numbers between 1 and N-1 which are relatively prime to N." << endl;
	cout << "Let's calculate phi(" << R << ")." << endl;
	cout << "Press Enter to Continue..." << endl;
	cin.ignore();

	int phiR = (U-1)*(V-1);

	cout << "phi(" << R << ") equals " << phiR << "." << endl;
	cout << endl;
	cout << "Now, we can find valid values for our public (P and R) and private (Q and R) keys." << endl;
	cout << "In order to be valid, P and Q need to be relatively prime to phi(R)." << endl;
	cout << "Additionally, P and Q must have this relationship with phi(R): (P*Q) = 1 (mod phi(R))." << endl;
	cout << "To do this, we'll iterate through prime values for P that are not factors of phi(R)." << endl;
	cout << "Then, we'll iterate through integers to find the modular multiplicative inverse for Q." << endl;
	cout << "If we can't find a modular multiplicative inverse for Q within a reasonable number of iterations, we'll continue to iterate P." << endl;	
	cout << "Press Enter to Continue..." << endl;
	cin.ignore();

	// Find primes of phiR
	vector<int> phiR_primes;
    phiR_primes.push_back(2);
    phiR_primes.push_back(3);
	primes(phiR_primes,phiR);
	vector<int> phiR_factors = factors(phiR_primes,phiR); // Find factors of phiR
	//for(int i=0; i < phiR_factors.size(); i++){
	//	cout << phiR_factors[i] << " ";
	//}

	// Calculate P and Q
	int P, Q;
	bool abort = false, advance = false;
	for (int P_index=0; P_index<phiR_primes.size(); P_index++) // Iterate through primes smaller than phi(R)
	{
		// Iterate through factors
		for (int factor_index=0; factor_index<phiR_factors.size(); factor_index++)
		{
			if(phiR_primes[P_index] == phiR_factors[factor_index])
			{
				//cout << "Prime, " << phiR_primes[P_index] << ", is equivalent to factor: " << phiR_factors[factor_index] << endl; 
				advance = true;
				break;
			}
		}

		if (advance) // Advance if prime is the same as a factor 
		{
			advance = false;
		} 
		else
		{
			// Iterate through K to find Q
			int K = 1;
			while (K < 500) // Limit to 500 iterations
			{
				/*
					Calculate Q, such that (P*Q) = 1 (mod phi(R))
					Then, find a value for P*Q that is congruent to: 1 (mod phi(R))
					Do this by multiplying phi(R) by a positive integer, K, then adding 1
					The resulting value, "temp", gives a value to solve Q = temp/P
					With a hypothetical value for P, phiR_primes[P_index], check if Q is an integer					
				*/
				int temp = ((K*phiR) + 1); 
				if (temp % phiR_primes[P_index] == 0) // If Q is an integer
				{
					Q = temp/phiR_primes[P_index];
					P = phiR_primes[P_index];
					abort = true;
					break;
				}
				else
				{
					K+=1;
				}
			}
		}

		if (abort)
		{
			break;
		}
	}

	cout << "Our value for P, is " << P << "." << endl;
	cout << "Our value for Q, is " << Q << "." << endl;
	cout << "Now, we need to generate an integer to character mapping, T, for each letter in the alphabet." << endl;
	cout << "However, values of T must always be less than R, relatively prime to R, and cannot be 1." << endl;
	cout << "Press Enter to Continue..." << endl;
	cin.ignore();

	// Generate vector for uppercase alphabet
	vector<char> alphabet;
	get_alphabet(alphabet);

	// Generate T values for each letter in the alphabet
	/*
		Assuming T will always be smaller than R: alphabet.size() < # of primes(U*V)
		This is why U and V must be 7 or larger.
	*/
	cout << "Our (26) values for T are... " << endl;
	vector<int> T;
	int T_iter = 1; // Will not be 1, Iterate to 2
	while(T.size() < alphabet.size())
	{
		T_iter++;
		if (check_rel_prime(T_iter,R)){ // Check if relatively prime
			cout << T_iter << " ";
			T.push_back(T_iter);
		}
	}
	cout << endl;

	// Get user word, convert with character-to-integer, T, mapping
	cout << "Now, we need a message." << endl;
	cout << "Please enter a word between 1 and 10 characters long: " << endl;
	string M = "";
	while (M.size() < 1 || M.size() > 10)
	{
		getline(cin,M);
		if(M.size() < 1 || M.size() > 10)
		{
			cout << "Enter a word between 1 and 10 characters long." << endl;
		} else {
			for(int i=0;i<M.size();i++)
			{
				if (!isalpha(M[i]))
				{
					cout << "Only use alphabetic letters." << endl;
					M = "";
					break;
				}
			}
		}
	}

	locale loc;
	vector<int> M_T;
	for (string::size_type i=0; i<M.length(); ++i){ // Iterate by letter
		M[i] = toupper(M[i],loc); // Convert to uppercase
		vector<char>::iterator alpha_p = find(alphabet.begin(),alphabet.end(),M[i]); // Find index in alphabet
		int index = alpha_p - alphabet.begin(); // From beginning
		M_T.push_back(T[index]); // Integer in T mapping
	}

	cout << "Your word is " << M << "." << endl;
	cout << "Converted into T values, your word is now: " << endl;
	for(int i=0; i < M_T.size(); i++){
		cout << M_T[i] << " ";
	}
	cout << endl;
	cout << endl;

	// Encrypt and Decrypt
	cout << "We can now encrypt the message with the public key values, P and R." << endl;
	cout << "The equation to do this is: X = (T^P) % R, for each value of T." << endl;
	cout << "Remember, P is " << P << " and R is " << R << "." << endl;
	cout << "Press Enter to Continue..." << endl;
	cin.ignore();

	cout << "Your encrypted message, X, is..." << endl;

	vector<mpz_class> X;
	for(int i=0; i < M_T.size(); i++){
		mpz_class p(P), t(M_T[i]), r(R), X_i; // Convert to multi-precision classes
		mpz_powm(X_i.get_mpz_t(), t.get_mpz_t(), p.get_mpz_t(), r.get_mpz_t()); // X = (T^P) % R
		cout << X_i << " ";
		X.push_back(X_i);
	}
	cout << endl;
	cout << endl;	
	cout << "Let's decrypt the message with our private keys, Q and R." << endl;
	cout << "The equation to do this is: T = (X^Q) % R, for each value of X." << endl;
	cout << "Remember, Q is " << Q << " and R is " << R << "." << endl;
	cout << "Press Enter to Continue..." << endl;
	cin.ignore();

	cout << "Your decrypted message, T, is..." << endl;

	vector<mpz_class> T2;
	for(int i=0; i < M_T.size(); i++){
		mpz_class q(Q), r(R), T2_i; // Convert to multi-precision classes
		mpz_powm(T2_i.get_mpz_t(), X[i].get_mpz_t(), q.get_mpz_t(), r.get_mpz_t()); // T = (X^Q) % R
		cout << T2_i << " ";
		T2.push_back(T2_i);
	}
	cout << endl;
	cout << endl;

	// Convert back to characters
	cout << "Let's convert this back into letters using the character mapping." << endl;
	cout << "Your original word was... ";
	for(int i=0; i < M_T.size(); i++){
		unsigned int T_i = T2[i].get_ui(); // T2[i].get_ui converts mpz_class to unsigned integer
		if (T_i > numeric_limits<unsigned int>::max()){ // Safety in case conversion from mpz_class produces a number too large
			break;
		}
		vector<int>::iterator T_p = find(T.begin(),T.end(),T_i); // Find index in T
		int index = T_p - T.begin();
		cout << alphabet[index];
	}
	cout << "!" << endl;
	
	return 0;
}