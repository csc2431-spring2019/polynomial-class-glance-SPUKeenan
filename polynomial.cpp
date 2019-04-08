#include "polynomial.h"

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cfloat>
#include <cmath>

using std::pow;
using std::istream;
using std::ostream;
using std::string;
using std::stringstream;
using std::fixed;
using std::setprecision;
using std::showpos;

// Produces: creates a dynamically allocated array
// Params: size_t degree, carries the degree size
// Returns: nothing
// Format Error: none
Polynomial::Polynomial(size_t degree) : _degree(degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = 0.0;
	}
}
// Produces: a deep copy of the array
// Params: size_t degree, degree of array, const float* coefficients, array
// Returns: nothing
// Format Error: none
Polynomial::Polynomial(size_t degree, const float* coefficients): _degree(degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = coefficients[i];
	}
}
// Produces: assigns values to the array
// Params: const Polynomial& polynomial, polynomial access variable
// Returns: nothing
// Format Error: none
Polynomial::Polynomial(const Polynomial& polynomial): _degree(polynomial._degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = polynomial._coefficients[i];
	}
}
// Produces: destructor for dynamically allocated memory
// Params: none
// Returns: nothing
// Format Error: none
Polynomial::~Polynomial(){
		delete []_coefficients;
}
// Produces: the sum of two polynomials
// Params: const Polu	& rhs, the second arrayt of data that we are adding
// Returns: the max degree and the solution
// Format Error: none
const Polynomial Polynomial::Sum(const Polynomial& rhs)const{
	size_t degreeMax = 0;
	size_t degreeMin = 0;

	if (_degree > rhs._degree){
		degreeMax = _degree;
		degreeMin = rhs._degree;
	} else{
		degreeMax = rhs._degree;
		degreeMin = _degree;
	}
	float solution[degreeMax + 1];
	if(_degree > rhs._degree){
		for (size_t i = 0; i < degreeMax + 1; i++) {
			solution[i] = _coefficients[i];
		}
		for (size_t i = 0; i < degreeMin + 1; i++) {
			solution[i] += rhs._coefficients[i];
		}
	}
	else {
		for (size_t i = 0; i < degreeMax; i++) {
			solution[i] = rhs._coefficients[i];
		}
		for (size_t i = 0; i < degreeMin + 1; i++) {
			solution[i] += _coefficients[i];
		}
	}
	return Polynomial(degreeMax, solution);
}
// Produces: subtracts two polynomials from eachother
// Params: const Polynomial& rhs, second array
// Returns: max degree and the solution array
// Format Error: none
const Polynomial Polynomial::Subtract(const Polynomial& rhs)const{
	size_t degreeMax = 0;
	size_t degreeMin = 0;

	if (_degree > rhs._degree){
		degreeMax = _degree;
		degreeMin = rhs._degree;
	} else{
		degreeMax = rhs._degree;
		degreeMin = _degree;
	}
	float solution[degreeMax + 1];
	if(_degree > rhs._degree){
		for (size_t i = 0; i < degreeMax + 1; i++) {
			solution[i] = _coefficients[i];
		}
		for (size_t i = 0; i < degreeMin + 1; i++) {
			solution[i] -= rhs._coefficients[i];
		}
	}
	else {
		for (size_t i = 0; i < degreeMax; i++) {
			solution[i] = rhs._coefficients[i];
		}
		for (size_t i = 0; i < degreeMin + 1; i++) {
			solution[i] -= _coefficients[i];
		}
	}
	return Polynomial(degreeMax, solution);
}
// Produces: a negative version of the polynomial array
// Params: none
// Returns: polynomial class retVal
// Format Error: none
const Polynomial Polynomial::Minus()const{
	Polynomial retVal(*this);
	for (size_t i = 0; i < _degree + 1; i++) {
		retVal._coefficients[i] *= -1;
	}
	return retVal;
}
// Produces: Multiplies two polynomials
// Params: const Polynomial& rhs, the second array to multiply with
// Returns: the max degree and the solution
// Format Error: none
const Polynomial Polynomial::Multiply(const Polynomial& rhs)const{
	size_t newDegree = _degree + rhs._degree;
	float solution[newDegree + 1];
	for (size_t i = 0; i < newDegree; i++) {
		solution[i] = 0.0;
	}
	for (size_t i = 0; i <= _degree; i++) {
		for (size_t j = 0; j <= rhs._degree; j++) {
			solution[i+j] += _coefficients[i] * rhs._coefficients[j];
		}
	}
	return Polynomial(newDegree, solution);
}
// Produces: Would Divides two arrays
// Params: const Polynomial& rhs, second array
// Returns: nothing
// Format Error: none
const Polynomial Polynomial::Divide(const Polynomial& rhs)const{
	return Polynomial(0);
}
// Produces: takes the derivitive of a polynomial
// Params: none
// Returns: the max degree and the solution
// Format Error: none
const Polynomial Polynomial::Derive()const{
	float solution [_degree];
	size_t newDegree = _degree - 1;

	for (size_t i = 0; i < _degree; i++) {
		solution[i] = 0.0;
	}
	for (size_t i = 0; i < _degree; i++) {
		solution[i] = (i+1) * _coefficients[i+1];
	}
	return Polynomial(newDegree, solution);
}
// Produces: evaluates the polynomial with a value of x
// Params: float x, the number we are evauluating at
// Returns: the solution
// Format Error: none
float Polynomial::Evaluate(float x)const{
	float solution [_degree + 1];
	float answer = 0.0;
	for (size_t i = 0; i < _degree + 1; i++){
		solution[i] = _coefficients[i];
	}
	for (size_t j = 0; j < _degree + 1; j++) {
		answer += solution[j] * (pow(x, j));
	}
	return answer;
}
// Produces: the integral of the polynomial
// Params: float start, the starting point, float end, the end point
// Returns: the definite integral
// Format Error: none
float Polynomial::Integrate(float start, float end)const{
	float solution[_degree + 2];
	float endAnswer = 0.0, startAnswer = 0.0, answer = 0.0;
	for (size_t i = 0; i < _degree + 1; i++) {
		solution[i] = _coefficients[i];
	}
	for (size_t j = 0; j < _degree + 1; j++) {
		endAnswer += ((solution[j] / (j + 1)) * pow(end, (j + 1)));
		startAnswer += ((solution[j] / (j + 1)) * pow(start, (j + 1)));
	}
	answer = endAnswer - startAnswer;
	return answer;
}
// Produces: allows us to use '=' as an operator
// Params: const Polynomial& rhs
// Returns: this, the polynomial
// Format Error: none
const Polynomial& Polynomial::operator=(const Polynomial& rhs){
	if (&rhs == this){
		return *this;
	}
	if (_degree != rhs._degree){
		if (_coefficients){
			delete[] _coefficients;
		}
		_degree = rhs._degree;
		_coefficients = new float[_degree + 1];
	}
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = rhs._coefficients[i];
	}
	return *this;
}
// Produces: checks to make sure that the answer we get and what it should be are actually equal
// Params: const Polynomial& rhs, second array
// Returns: true or false
// Format Error: none
bool Polynomial::Equals(const Polynomial& rhs)const{
	if (_degree != rhs._degree){
		return false;
	}
	for (size_t i=0; i < _degree; i++){
		if (abs(_coefficients[i] - rhs._coefficients[i]) > 0.0001){
			return false;
		}
	}
	return true;
}
// Produces: an output in string format
// Params: none
// Returns: the string we are outputting
// Format Error: none
string Polynomial::ToString()const{
	stringstream ss;
	for (size_t i = _degree; i > 0; i--) {
		ss << showpos << fixed << setprecision(2) << _coefficients[i] << "x^" << i << " ";
	}
	ss << showpos << fixed << setprecision(2) << _coefficients[0];
	return ss.str();
}
// Produces: output of the data
// Params: ostream& output, output variable
// Returns: the output of the data
// Format Error: none
ostream& Polynomial::Write(ostream& output)const{
	output << _degree << " ";
	for (size_t i = 0; i < _degree + 1; i++) {
		output << _coefficients[i] << " ";
	}
	return output;
}
// Produces: reads data from the input stream
// Params: istream& input, input variable
// Returns: the input
// Format Error: none
istream& Polynomial::Read(istream& input){
	size_t degree;
	input >> degree;
	if (input.fail()){
		return input;
	}
	float* coefficients = new float[degree + 1];
	for (size_t i = 0; i < degree + 1; i++) {
		input >> coefficients[i];
		if (input.fail()){
			delete[] coefficients;
			return input;
		}
	}

	if (degree != _degree){
		if (_coefficients){
			delete[] _coefficients;
		}
		_degree = degree;
		_coefficients = coefficients;
	}else{
		for (size_t i = 0; i < _degree + 1; i++) {
			_coefficients[i] = coefficients[i];
		}
		delete[] coefficients;
	}
	return input;
}
