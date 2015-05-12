/*
 * BaseBox.cpp
 *
 *  Created on: 17/mar/2015
 *      Author: lorenzo
 */

#include "BaseBox.h"

template<typename number>
BaseBox<number>::BaseBox() {

}

template<typename number>
BaseBox<number>::~BaseBox() {

}

template class BaseBox<float>;
template class BaseBox<double>;
