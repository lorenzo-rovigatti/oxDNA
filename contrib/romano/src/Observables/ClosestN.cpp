/*
 * ClosestMap.cpp
 *
 * Created on Mar 25, 2019
 *       Author: Poppleton
 */

#include "ClosestN.h"
#include <sstream>

ClosestN::ClosestN() {
	_rmax = 5.;
	_nmax = 6;
}

ClosestN::~ClosestN() {

}

void ClosestN::init() {
	BaseObservable::init();
}

void ClosestN::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "nmax", &_nmax, 0);
}

neighs::neighs(int max_n) {
	_max_size = max_n;
	ds = std::vector<number>(_max_size, 1.e9);
	ns = std::vector<int>(_max_size, -1);
}

neighs::~neighs() {

}

void neighs::add_pair(int i, number d) {
	//std::cout << "entering; " << i << " " << d << std::endl;
	std::vector<number>::iterator itd = ds.begin();
	std::vector<int>::iterator itn = ns.begin();
	while (itd != ds.end() && *itd < d) {
		itd ++;
		itn ++;
	}
	//for (auto e : ds) std::cout << e << " ";
	//for (auto e : ns) std::cout << e << " ";
	//std::cout << std::endl;
	ds.insert(itd,d);
	ns.insert(itn,i);
	ds.resize(_max_size);
	ns.resize(_max_size);
	//for (auto e : ds) std::cout << e << " ";
	//for (auto e : ns) std::cout << e << " ";
	//std::cout << std::endl;
}


std::string ClosestN::get_output_string(llint curr_step) {
	int N = _config_info->N();

	std::vector<neighs *> close;
	close.resize(N);
	for(int i = 0; i < N; i++) {
		close[i] = new neighs(_nmax);
	}

	for(int i = 0; i < N; i++) {
		for(int j = i + 1; j < N; j++) {
			LR_vector dist;
			LR_vector p1_com, p2_com;
			p1_com = _config_info->box->get_abs_pos(_config_info->particles()[i]);
			p2_com = _config_info->box->get_abs_pos(_config_info->particles()[j]);
			
			number myd2 = _config_info->box->min_image(p1_com, p2_com).norm();

			close[i]->add_pair(j, myd2);
			close[j]->add_pair(i, myd2);

			//if (i == 0) std::cin.get();
		}
	}

	std::stringstream outstr;
	int n = 0;
	for (auto c : close) {
		outstr << n << " - ";
		for (auto e : c->ns) {
			outstr << " " << e;
		}
		outstr << std::endl;
		n++;
	}

	return outstr.str();
}

/*
std::string ClosestN::get_output_string(llint curr_step) {
	int N = _config_info->N();

	std::vector<std::vector<int> > maxn;
	std::vector<std::vector<number> > close;

	maxn.resize(N);
	close.resize(N);
	for(int i = 0; i < N; i++) {
		close[i].resize(0);
		maxn[i].resize(0);
	}

	for(int i = 0; i < N; i++) {
		for(int j = i + 1; j < N; j++) {
			LR_vector dist;
			LR_vector p1_com, p2_com;
			p1_com = _config_info->box->get_abs_pos(_config_info->particles()[i]);
			p2_com = _config_info->box->get_abs_pos(_config_info->particles()[j]);
			
			number myd2 = _config_info->box->min_image(p1_com, p2_com).norm();

			close[i].push_back(myd2);
			maxn[i].push_back(j);
			if ((int)close[i].size() > _nmax) {
				// find largest

				//std::cout << "BEFORE " << std::endl;
				//for (int k = 1; k < _nmax + 1; k ++) {
				//	std::cout << maxn[i][k] << " " << close[i][k] << "  ";
				//}
				//std::cout << std::endl;
				number mymax = close[i][0];
				int nmymax = 0;
				for (int k = 1; k <= _nmax; k ++) {
					if (close[i][k] > mymax) {
						mymax = close[i][k];
						nmymax = k;
					}
				}
				//std::cout << "removing " << nmymax << " " << maxn[i][nmymax] << " " << close[i][nmymax] << std::endl;
				close[i].erase(close[i].begin()+nmymax);
				maxn[i].erase(maxn[i].begin()+nmymax);
				//std::cout << "AFTER " << std::endl;
				//for (int k = 1; k < _nmax; k ++) {
				//	std::cout << maxn[i][k] << " " << close[i][k] << "  ";
				//}
				//std::cout << std::endl;
				//std::cout << std::endl;
			}
			
			close[j].push_back(myd2);
			maxn[j].push_back(i);
			if ((int)close[j].size() > _nmax) {
				number mymax = close[j][0];
				int nmymax = 0;
				for (int k = 1; k <= _nmax; k ++) {
					if (close[j][k] > mymax) {
						mymax = close[j][k];
						nmymax = k;
					}
				}
				close[j].erase(close[j].begin()+nmymax);
				maxn[j].erase(maxn[j].begin()+nmymax);
			}
		}
	}

	//for(int i = 0; i < N; i++) std::cout << maxn[i].size(); 
	
	std::stringstream outstr;
	for(int i = 0; i < N; i++) {
		outstr << i << " -";
		for (auto j : maxn[i])
			outstr << " " << j;
		outstr << std::endl;
	}

	return outstr.str();
}
*/
