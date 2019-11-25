// GNU General Public License v3.0
//
// Author: Nicolas Salmieri Nov 2019

#pragma once

#include <vector>
#include <math.h>
#include <algorithm>
#include <stdint.h>
#include <assert.h>
#include <iostream>
#include <iomanip>

struct Np{
	template<typename T>
	struct NdArray:public std::vector<T>{

		NdArray():
			std::vector<T>(),dimentions({}){}
		NdArray(const NdArray& other):
			std::vector<T>(*&other),dimentions(other.dimentions){}

		NdArray operator*(const NdArray& other) const{
			assert(other.dimentions==dimentions);
			NdArray ret(*this);
			for (size_t i = 0; i < this->size(); ++i) {
				ret[i]*=other.at(i);
			}
			return ret;
		}

		NdArray pow(double power){
			NdArray ret(*this);
			for (auto& elem : ret) {
				elem=std::pow(elem,power);
			}
			return ret;
		}

		NdArray flatten(char order='C'){
			auto unused=order;unused++;
			NdArray ret(*this);
			ret.dimentions.resize(1,ret.size());
			return ret;
		}

		T& operator[](std::vector<int> slice){
			int offset=slice[0];
			for (size_t dim = 1; dim < slice.size(); ++dim) {
				offset+=dimentions[dim-1]*slice[dim];
			}
			return (*this)[offset];
		}
		T& operator[](int index){
			return std::vector<T>::operator[](index);
		}
		std::vector<size_t> dimentions=std::vector<size_t>();//! dimentions of the array
	};



	//! return an array filled of ones of dimentions dimentions
	NdArray<double> ones(std::vector<size_t> dimentions){
		int size=1;
		for (auto& dim : dimentions) {
			size*=dim;
		}

		NdArray<double> ret;
		ret.resize(size,1);
		ret.dimentions=dimentions;

		return ret;
	}
	//! return an array filled of zeros of dimentions dimentions
	NdArray<double> zeros(std::vector<size_t> dimentions){
		int size=1;
		for (auto& dim : dimentions) {
			size*=dim;
		}

		NdArray<double> ret;
		ret.resize(size,0);
		ret.dimentions=dimentions;

		return ret;
	}
	NdArray<double> arange(double start,double stop,double step=1){
		NdArray<double> ret;
		ret.reserve((stop-start)/step+2/*+2 just in case*/);

		for (double i = start; i < stop; i+=step) {
			ret.push_back(i);
		}
		ret.pop_back();

		ret.dimentions.reserve(1);
		ret.dimentions.resize(1,ret.size());
		return ret;
	}
	std::vector<size_t> size(const NdArray<double>& array){
		return array.dimentions;
	}
	NdArray<double> exp(const NdArray<double>& array){
		NdArray<double> ret(array);
		for (auto& elem : ret) {
			elem=std::exp(elem);
		}
		return ret;
	}
	double power(double v,double power){
		return std::pow(v,power);
	}
	double pi=3.14159265359;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const Np::NdArray<T>& array){
	for (auto& elem : array) {
		os << elem << ' ';
	}

	return os;
}

template<typename T>
Np::NdArray<T> operator*(double factor,const Np::NdArray<T>& array){
	auto ret=Np::NdArray<T>(array);
	for (auto& elem : ret) {
		elem*=factor;
	}
	return ret;
}
template<typename T>
Np::NdArray<T> operator*(const Np::NdArray<T>& array,double factor){
	return factor*array;
}

template<typename T>
Np::NdArray<T> operator-(double value,const Np::NdArray<T>& array){
	auto ret=Np::NdArray<T>(array);
	for (auto& elem : ret) {
		elem=value-elem;
	}
	return ret;
}
template<typename T>
Np::NdArray<T> operator-(const Np::NdArray<T>& array,double value){
	auto ret=Np::NdArray<T>(array);
	for (auto& elem : ret) {
		elem-=value;
	}
	return ret;
}

template<typename T>
T min(const Np::NdArray<T>& array){
	return *std::min_element(array.begin(),array.end());
}

template<typename T>
T max(const Np::NdArray<T>& array){
	return *std::max_element(array.begin(),array.end());
}

template<typename E>
void print(E v) {
	std::cout << v << std::endl;
}
template<typename E, typename... Args>
void print(E first, Args... args) {
	std::cout << first <<" ";
	print(args...);
}

std::vector<int> range(int start,int stop,int step=1){
	std::vector<int> ret;
	for (int i = start; i < stop; i+=step) {
		ret.push_back(i);
	}
	return ret;
}



