/*
* contact_map.cpp
*
* Created on Mar 25, 2019
*       Author: Poppleton
*/

#include "ContactMap.h"
#include <sstream>

template<typename number>
ContactMap<number>::ContactMap() {

}

template<typename number>
ContactMap<number>::~ContactMap() {

}

template<typename number>
void ContactMap<number>::init(ConfigInfo<number> &config_info) {
    BaseObservable<number>::init(config_info);
}

template<typename number>
std::string ContactMap<number>::get_output_string(llint curr_step) {
    int n = *this -> _config_info.N;
    int s = ((n*n)-n)/2;
    number *cmap  = new number[s];
    
    int k=0;
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            LR_vector<number> dist;
            LR_vector<number> p1_com, p2_com;
            p1_com = this->_config_info.box->get_abs_pos(this -> _config_info.particles[i]);
            p2_com = this->_config_info.box->get_abs_pos(this -> _config_info.particles[j]);
            cmap[k] = this->_config_info.box->min_image(p1_com,p2_com).module();
            k++;
        }
    }
    std::stringstream outstr;
    for (int i =0; i < s; i++) {
        outstr << cmap[i];
        outstr << ' ';
    }

    delete [] cmap;
    return outstr.str();
}


template class ContactMap<float>;
template class ContactMap<double>;
