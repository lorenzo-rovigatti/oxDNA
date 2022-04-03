#include "GaussTrapMetaRatio1D.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

#include <string>
#include <sstream>
#include <iostream>

std::vector<int> split3(const std::string& s, char delimiter){
   std::vector<int> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(std::stoi(token)); 
   }
   return tokens;
}

std::vector<double> split3_potential_line(const std::string& s){
   std::vector<double> potential_line;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, ','))
   {
      potential_line.push_back(std::stod(token)); 
   }
   return potential_line;
}

inline double interpolatePotential2(double a,double dX,double xmin,const std::vector<double > & potential_grid){
    
    double a_left  = dX*std::floor(a/dX);
    double a_right = a_left + dX;

    double ia_left  = std::floor((a-xmin)/dX);
    double ia_right = ia_left + 1;

    double f1 = potential_grid[ia_left];
    double f2 = potential_grid[ia_right];

    double f =  (a_right - a) /dX * f1 + (a - a_left) /dX * f2; 
    return f;
}

inline double get_x_force2(const double x,const double y, const double dX,const double xmin,const std::vector<double> & potential_grid){
    const double delta = dX/5.;
    const double ratio = y / x;
    const double a = std::atan(ratio);
    const double Uprime = (interpolatePotential2(a+delta,dX,xmin,potential_grid) - interpolatePotential2(a-delta,dX,xmin,potential_grid) ) / (2*delta);
    const double prefactor = - (1./(1+std::pow(ratio,2)))*(y/std::pow(x,2)); // negative sign comes from derivative wrt denominator
    const double force = - prefactor * Uprime;
    return force;
}
inline double get_y_force2(const double x,const double y, const double dX,const double xmin,const std::vector<double> & potential_grid){
    const double delta = dX/5.;
    const double ratio = y / x;
    const double a = std::atan(ratio);
    const double Uprime = (interpolatePotential2(a+delta,dX,xmin,potential_grid) - interpolatePotential2(a-delta,dX,xmin,potential_grid) ) / (2*delta);
    const double prefactor = (1./(1+std::pow(ratio,2)))*(1/x);
    const double force = - prefactor * Uprime;
    return force;
}


template<typename number>
GaussTrapMetaRatio1D<number>::GaussTrapMetaRatio1D() : BaseForce<number>() {
	_ref_id = -2;
    _p1a = {};
    _p2a = {};
	_p1a_ptr = {};
	_p2a_ptr = {};
    _p1b = {};
    _p2b = {};
	_p1b_ptr = {};
	_p2b_ptr = {};

    xmin = 0;
    xmax = 10;
    dX = 0.1;

    PBC = false;
	_box_ptr = NULL;
    _mode = -1;
}

template <typename number>
void GaussTrapMetaRatio1D<number>::get_settings (input_file &inp) {

    std::string _p1a_string;
    std::string _p2a_string;
    std::string _p1b_string;
    std::string _p2b_string;

	getInputString (&inp, "p1a", _p1a_string, 1);
	getInputString (&inp, "p2a", _p2a_string, 1);
	getInputString (&inp, "p1b", _p1b_string, 1);
	getInputString (&inp, "p2b", _p2b_string, 1);

    _p1a = split3(_p1a_string,',');
    _p2a = split3(_p2a_string,',');
    _p1b = split3(_p1b_string,',');
    _p2b = split3(_p2b_string,',');

	getInputInt (&inp, "mode", &this->_mode, 1);

	getInputBool (&inp, "PBC", &PBC, 0);

	getInputNumber(&inp, "xmin", &this->xmin, 1);
	getInputNumber(&inp, "xmax", &this->xmax, 1); // these are both inclusive
	getInputInt(&inp, "N_grid", &this->N_grid, 1); // we don't want to take this directly.

    std::string potential_string;
	getInputString (&inp, "potential_grid", potential_string, 1);
    potential_grid = split3_potential_line(potential_string); 
}

template <typename number>
void GaussTrapMetaRatio1D<number>::init (BaseParticle<number> ** particles, int N, BaseBox<number> * box_ptr){
	this->_box_ptr = box_ptr;

    for (int i = 0; i < _p1a.size(); i++){
        if (_p1a[i] < 0 || _p1a[i] >= N) throw oxDNAException ("Invalid reference particle %d for Gauss Trap", _p1a[i]);
    }
    for (int i = 0; i < _p2a.size(); i++){
        if (_p2a[i] < 0 || _p2a[i] >= N) throw oxDNAException ("Invalid reference particle %d for Gauss Trap", _p2a[i]);
    }

    for (int i = 0; i < _p1a.size(); i++){
        _p1a_ptr.push_back(particles[_p1a[i]]);
    }
    for (int i = 0; i < _p2a.size(); i++){
        _p2a_ptr.push_back(particles[_p2a[i]]);
    }
    for (int i = 0; i < _p1b.size(); i++){
        _p1b_ptr.push_back(particles[_p1b[i]]);
    }
    for (int i = 0; i < _p2b.size(); i++){
        _p2b_ptr.push_back(particles[_p2b[i]]);
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
        for (int i = 0; i < _p1b.size(); i++){
            particles[_p1b[i]]->add_ext_force(this); 
        }
    }
    else if (this->_mode == 4){
        for (int i = 0; i < _p2b.size(); i++){
            particles[_p2b[i]]->add_ext_force(this); 
        }
    }
    this->dX = (this->xmax - this->xmin) / (this->N_grid-1);
}

template<typename number>
LR_vector<number> GaussTrapMetaRatio1D<number>::_distance(LR_vector<number> u, LR_vector<number> v) {
	if (this->PBC) return this->_box_ptr->min_image(u, v);
	else return v - u;
}

template<typename number>
LR_vector<number> GaussTrapMetaRatio1D<number>::value (llint step, LR_vector<number> &pos) {

    LR_vector<number> p1a_vec = {0,0,0};
    LR_vector<number> p2a_vec = {0,0,0};
    LR_vector<number> p1b_vec = {0,0,0};
    LR_vector<number> p2b_vec = {0,0,0};

    for (int i = 0 ; i < _p1a_ptr.size(); i++){
        p1a_vec += this->_box_ptr->get_abs_pos(_p1a_ptr[i]);  
    }
    p1a_vec = p1a_vec / (double)_p1a_ptr.size();

    for (int i = 0 ; i < _p2a_ptr.size(); i++){
        p2a_vec += this->_box_ptr->get_abs_pos(_p2a_ptr[i]);  
    }
    p2a_vec = p2a_vec / (double)_p2a_ptr.size();

    for (int i = 0 ; i < _p1b_ptr.size(); i++){
        p1b_vec += this->_box_ptr->get_abs_pos(_p1b_ptr[i]);  
    }
    p1b_vec = p1b_vec / (double)_p1b_ptr.size();

    for (int i = 0 ; i < _p2b_ptr.size(); i++){
        p2b_vec += this->_box_ptr->get_abs_pos(_p2b_ptr[i]);  
    }
    p2b_vec = p2b_vec / (double)_p2b_ptr.size();

	LR_vector<number> dra = this->_distance(p2a_vec, p1a_vec);
	LR_vector<number> drb = this->_distance(p2b_vec, p1b_vec);

    double my_x = dra.module();
    double my_y = drb.module();

    double angle = std::atan(my_y/my_x);
    int ix_left = std::floor((angle-this->xmin)/this->dX);
    int ix_right = ix_left+1;

    double meta_Fx = 0;
    double meta_Fy = 0;

    if ((ix_left < 0) || (ix_right > N_grid-1) ){
        std::cout << "off grid!" << std::endl;
    }

    else{
        meta_Fx = get_x_force2(my_x,my_y,dX,xmin,potential_grid);
        meta_Fy = get_y_force2(my_x,my_y,dX,xmin,potential_grid);
        //std::cout << meta_Fx << std::endl;
        //std::cout << meta_Fy << std::endl;
    }

    const LR_vector<number> accumulated_force_x = dra*(meta_Fx / dra.module()) ;
    const LR_vector<number> accumulated_force_y = drb*(meta_Fy / drb.module()) ;

    // deliberately wrong sign for experiment
    if (this->_mode == 1){
        return accumulated_force_x / _p1a.size(); 
    }
    else if (this->_mode == 2){
        return (accumulated_force_x*-1) / _p2a.size();
    }
    // could one of these signs be wrong?
    else if (this->_mode == 3){
        return accumulated_force_y / _p1b.size();  
    }
    else if (this->_mode == 4){
        return (accumulated_force_y*-1) / _p2b.size();
    }
}


template<typename number>
number GaussTrapMetaRatio1D<number>::potential (llint step, LR_vector<number> &pos) {



    LR_vector<number> p1a_vec = {0,0,0};
    LR_vector<number> p2a_vec = {0,0,0};
    LR_vector<number> p1b_vec = {0,0,0};
    LR_vector<number> p2b_vec = {0,0,0};

    for (int i = 0 ; i < _p1a_ptr.size(); i++){
        p1a_vec += this->_box_ptr->get_abs_pos(_p1a_ptr[i]);  
    }
    p1a_vec = p1a_vec / (double)_p1a_ptr.size();

    for (int i = 0 ; i < _p2a_ptr.size(); i++){
        p2a_vec += this->_box_ptr->get_abs_pos(_p2a_ptr[i]);  
    }
    p2a_vec = p2a_vec / (double)_p2a_ptr.size();

    for (int i = 0 ; i < _p1b_ptr.size(); i++){
        p1b_vec += this->_box_ptr->get_abs_pos(_p1b_ptr[i]);  
    }
    p1b_vec = p1b_vec / (double)_p1b_ptr.size();

    for (int i = 0 ; i < _p2b_ptr.size(); i++){
        p2b_vec += this->_box_ptr->get_abs_pos(_p2b_ptr[i]);  
    }
    p2b_vec = p2b_vec / (double)_p2b_ptr.size();

	LR_vector<number> dra = this->_distance(p2a_vec, p1a_vec);
	LR_vector<number> drb = this->_distance(p2b_vec, p1b_vec);

    double my_x = dra.module();
    double my_y = drb.module();
    double angle = my_y/my_x;
    int ix_left = std::floor((angle-this->xmin)/this->dX);
    int ix_right = ix_left+1;

    double meta_Fx = 0;
    double meta_Fy = 0;

    if ((ix_left < 0) || (ix_right > N_grid-1) ){
        std::cout << "off grid!" << std::endl;
    }

    double my_potential = 0;

    const double ratio = my_y / my_x;
    const double a = std::atan(ratio);

    my_potential = interpolatePotential2(a,dX,xmin,potential_grid);

    int total_factor = _p1a_ptr.size() + _p2a_ptr.size() + _p1b_ptr.size() + _p2b_ptr.size();

	return my_potential / (number)(total_factor);
}

template class GaussTrapMetaRatio1D<double>;
template class GaussTrapMetaRatio1D<float>;











