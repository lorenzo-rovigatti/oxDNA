#include "GaussTrapAngle.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

#include <string>
#include <sstream>
#include <iostream>

std::vector<int> split4(const std::string& s, char delimiter){
   std::vector<int> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(std::stoi(token)); 
   }
   return tokens;
}

std::vector<double> split4_potential_grid(const std::string& s){
   std::vector<double> potential_grid;
   std::string token;
   std::istringstream tokenStream(s);

   while (std::getline(tokenStream, token, ','))
   {
      potential_grid.push_back(std::stod(token)); 
   }
   return potential_grid;
}

inline double interpolatePotential(const double my_x, const double dX, const double xmin,const std::vector<double> & potential_grid){
    double x_left  = dX*std::floor(my_x/dX);
    double x_right = x_left + dX;
    double ix_left  = std::floor((my_x-xmin)/dX);
    double ix_right = ix_left + 1;
    double f11 = potential_grid[ix_left];
    double f21 = potential_grid[ix_right];
    double fx =  (x_right - my_x) /dX * f11 + (my_x - x_left) /dX * f21;  
    return fx;
}

inline double get_x_force(const double x,const double dX,const double xmin,const std::vector<double> & potential_grid){
    double ix_left  = std::floor((x-xmin)/dX);
    double ix_right = ix_left + 1;
    return -(potential_grid[ix_right] - potential_grid[ix_left])/dX;
}

template<typename number>
GaussTrapAngle<number>::GaussTrapAngle() : BaseForce<number>() {
	_ref_id = -2;
    _p1a = {};
    _p2a = {};
    _p3a = {};
	_p1a_ptr = {};
	_p2a_ptr = {};
	_p3a_ptr = {};

    xmin = 0;
    xmax = 10;
    dX = 0.1;

    PBC = false;
	_box_ptr = NULL;
    _mode = -1;
}

template <typename number>
void GaussTrapAngle<number>::get_settings (input_file &inp) {

    std::string _p1a_string;
    std::string _p2a_string;
    std::string _p3a_string;

	getInputString (&inp, "p1a", _p1a_string, 1);
	getInputString (&inp, "p2a", _p2a_string, 1);
	getInputString (&inp, "p3a", _p3a_string, 1);

    _p1a = split4(_p1a_string,',');
    _p2a = split4(_p2a_string,',');
    _p3a = split4(_p3a_string,',');

	getInputInt (&inp, "mode", &this->_mode, 1);

	getInputBool (&inp, "PBC", &PBC, 0);

	getInputNumber(&inp, "xmin", &this->xmin, 1);
	getInputNumber(&inp, "xmax", &this->xmax, 1); // these are both inclusive
	getInputInt(&inp, "N_grid", &this->N_grid, 1);

    std::string potential_string;
	getInputString (&inp, "potential_grid", potential_string, 1);
    potential_grid = split4_potential_grid(potential_string); 

}

template <typename number>
void GaussTrapAngle<number>::init (BaseParticle<number> ** particles, int N, BaseBox<number> * box_ptr){
	this->_box_ptr = box_ptr;

    for (int i = 0; i < _p1a.size(); i++){
        if (_p1a[i] < 0 || _p1a[i] >= N) throw oxDNAException ("Invalid reference particle %d for Gauss Trap", _p1a[i]);
    }
    for (int i = 0; i < _p2a.size(); i++){
        if (_p2a[i] < 0 || _p2a[i] >= N) throw oxDNAException ("Invalid reference particle %d for Gauss Trap", _p2a[i]);
    }
    for (int i = 0; i < _p3a.size(); i++){
        if (_p3a[i] < 0 || _p3a[i] >= N) throw oxDNAException ("Invalid reference particle %d for Gauss Trap", _p3a[i]);
    }

    for (int i = 0; i < _p1a.size(); i++){
        _p1a_ptr.push_back(particles[_p1a[i]]);
    }
    for (int i = 0; i < _p2a.size(); i++){
        _p2a_ptr.push_back(particles[_p2a[i]]);
    }
    for (int i = 0; i < _p3a.size(); i++){
        _p3a_ptr.push_back(particles[_p3a[i]]);
    }

    if (this->_mode == 1){
        for (int i = 0; i < _p1a.size(); i++){
            particles[_p1a[i]]->add_ext_force(this); 
        }
    }
    else if (this->_mode == 2){
        for (int i = 0; i < _p2a.size(); i++){
            particles[_p2a[i]]->add_ext_force(this); 
        }
    }
    else if (this->_mode == 3){
        for (int i = 0; i < _p3a.size(); i++){
            particles[_p3a[i]]->add_ext_force(this); 
        }
    }

    this->dX = (this->xmax - this->xmin) / (this->N_grid-1);
}

template<typename number>
LR_vector<number> GaussTrapAngle<number>::_distance(LR_vector<number> u, LR_vector<number> v) {
	if (this->PBC) return this->_box_ptr->min_image(u, v);
	else return v - u;
}

template<typename number>
LR_vector<number> GaussTrapAngle<number>::value (llint step, LR_vector<number> &pos) {

    LR_vector<number> p1a_vec = {0,0,0};
    LR_vector<number> p2a_vec = {0,0,0};
    LR_vector<number> p3a_vec = {0,0,0};

    for (int i = 0 ; i < _p1a_ptr.size(); i++){
        p1a_vec += this->_box_ptr->get_abs_pos(_p1a_ptr[i]);  
    }
    p1a_vec = p1a_vec / (double)_p1a_ptr.size();

    for (int i = 0 ; i < _p2a_ptr.size(); i++){
        p2a_vec += this->_box_ptr->get_abs_pos(_p2a_ptr[i]);  
    }
    p2a_vec = p2a_vec / (double)_p2a_ptr.size();

    for (int i = 0 ; i < _p3a_ptr.size(); i++){
        p3a_vec += this->_box_ptr->get_abs_pos(_p3a_ptr[i]);  
    }
    p3a_vec = p3a_vec / (double)_p3a_ptr.size();

	LR_vector<number> dra1 = this->_distance(p2a_vec, p1a_vec);
	LR_vector<number> dra2 = this->_distance(p2a_vec, p3a_vec);

    LR_vector<number> dra1_normed = dra1 / dra1.module();
    LR_vector<number> dra2_normed = dra2 / dra2.module();
    
    double dot_product = dra1_normed * dra2_normed;
    double angle = std::acos(dot_product);
    int ix_left = std::floor((angle-this->xmin)/this->dX);
    int ix_right = ix_left+1;

    double xforce = 0;

    if ((ix_left < 0) || (ix_right > N_grid-1) ){
        std::cout << "off grid!" << std::endl;
    }

    else{
        xforce = get_x_force(angle,dX,xmin,potential_grid);
    }

    // why is this necessary?
    double prefactor = - xforce/std::pow(1-dot_product,0.5);

    // this isn't fucking working

    if (this->_mode == 1){
        double r1_factor = - dot_product / dra1.module();
        double r2_factor =  1 / dra1.module(); 
        return ((dra1_normed*r1_factor) + (dra2_normed*r2_factor))*prefactor/(double)_p1a_ptr.size();
    }
    if (this->_mode == 3){
        double r2_factor = - dot_product / dra2.module();
        double r1_factor =  1 / dra2.module(); 
        return ((dra1_normed*r1_factor) + (dra2_normed*r2_factor))*prefactor/(double)_p3a_ptr.size();
    }
    if (this->_mode == 2){
        double r1_factor_A = - dot_product / dra1.module();
        double r2_factor_A =  1 / dra1.module(); 
        double r2_factor_B = - dot_product / dra2.module();
        double r1_factor_B =  1 / dra2.module(); 
        return ((dra1_normed*(-r1_factor_A-r1_factor_B)) + (dra2_normed*(-r2_factor_A-r2_factor_B)))*prefactor/(double)_p2a_ptr.size();
    }

}


template<typename number>
number GaussTrapAngle<number>::potential (llint step, LR_vector<number> &pos) {


    LR_vector<number> p1a_vec = {0,0,0};
    LR_vector<number> p2a_vec = {0,0,0};
    LR_vector<number> p3a_vec = {0,0,0};

    for (int i = 0 ; i < _p1a_ptr.size(); i++){
        p1a_vec += this->_box_ptr->get_abs_pos(_p1a_ptr[i]);  
    }
    p1a_vec = p1a_vec / (double)_p1a_ptr.size();

    for (int i = 0 ; i < _p2a_ptr.size(); i++){
        p2a_vec += this->_box_ptr->get_abs_pos(_p2a_ptr[i]);  
    }
    p2a_vec = p2a_vec / (double)_p2a_ptr.size();

    for (int i = 0 ; i < _p3a_ptr.size(); i++){
        p3a_vec += this->_box_ptr->get_abs_pos(_p3a_ptr[i]);  
    }
    p3a_vec = p3a_vec / (double)_p3a_ptr.size();

	LR_vector<number> dra1 = this->_distance(p2a_vec, p1a_vec);
	LR_vector<number> dra2 = this->_distance(p2a_vec, p3a_vec);

    LR_vector<number> dra1_normed = dra1 / dra1.module();
    LR_vector<number> dra2_normed = dra2 / dra2.module();
    
    double dot_product = dra1_normed * dra2_normed;
    double angle = std::acos(dot_product);
    int ix_left = std::floor((angle-this->xmin)/this->dX);
    int ix_right = ix_left+1;

    double my_potential = 0;

    if ((ix_left < 0) || (ix_right > N_grid-1) ){
        std::cout << "off grid!" << std::endl;
    }

    else{
        my_potential = interpolatePotential(angle, dX, xmin, potential_grid);
    }

    int total_factor = _p1a_ptr.size() + _p2a_ptr.size() + _p3a_ptr.size() ;

	return my_potential / (number)(total_factor);
}

template class GaussTrapAngle<double>;
template class GaussTrapAngle<float>;











