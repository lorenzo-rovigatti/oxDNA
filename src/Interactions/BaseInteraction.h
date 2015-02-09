/*
 * BaseInteraction is a class that is meant "contain" all the interaction
 * types
 */

#ifndef BASE_INTERACTION_H
#define BASE_INTERACTION_H

#include <map>
#include <cfloat>
#include <fstream>
#include <set>
#include <vector>

#include "../defs.h"
#include "../Particles/BaseParticle.h"
#include "../Utilities/Utils.h"
#include "../Utilities/oxDNAException.h"
#include "./Mesh.h"

using namespace std;

/**
 * @brief Abstract class defining the interaction interface.
 *
 * We implement the Curiously Recurring Template Pattern (see
 * http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) and therefore
 * we have to keep this class and BaseInteraction separated.
 */
template <typename number>
class IBaseInteraction {
protected:
	number _box_side;

	/// Is used by _create_cells() to avoid useless memory allocations.
	number _last_box_side;

	/// Keeps track of all the cells containing at least one particle.
	std::set<int> _not_empty_cells;

	/// This is useful for "hard" potentials
	bool _is_infinite;

	/// This is useful to create initial configurations. it is used by the generator functions
	number _energy_threshold;

	number _rcut, _sqr_rcut;

	char _topology_filename[256];

	int _cells_N_side;
	int *_cells_head;
	int *_cells_index;
	int *_cells_next;
	int **_cells_neigh;

	/**
	 * @brief Allocates and initializes cells used to compute interaction energies.
	 *
	 * The caller must take care of deleting the _cells_neigh array. The implementation of this
	 * method is complicated by the fact that it has to be O(N) and not O(box_side^3) for all
	 * box_side values. In order to have good performances we resort to bookkeeping. This is why
	 * we use _last_box_side and _not_empty_cells. If the last optional argument is true,
	 * the bookkeeping is not used and all the cells gets initialized. This brings back the
	 * O(box_side^3) computational complexity.
	 *
	 * @param particles
	 * @param N
	 * @param box_side
	 * @param init_all_neighs initializes all the cells. If true, then this method becomes O(box_side^3)
	 */
	void _create_cells(BaseParticle<number> **particles, int N, number box_side, bool init_all_neighs=false);

	/**
	 * @brief Delete all the cell-related arrays but _cells_neigh.
	 */
	void _delete_cells();

	/**
	 * @brief Free _cells_neigh.
	 *
	 * It uses _not_empty_cells to check whether _cells_neigh[i] is allocated or not. This method
	 * has to be called by all the methods which also call _create_cells to avoid memory leaks.
	 */
	void _delete_cell_neighs();
public:
	IBaseInteraction();
	virtual ~IBaseInteraction();

	virtual void set_box_side(number box_side) { _box_side = box_side; }

	virtual void get_settings(input_file &inp);

	/**
	 * Initialization of class constants.
	 */
	virtual void init() = 0;

	/**
	 * @brief Handles particle allocation. Child classes must implement it.
	 *
	 * @param particles
	 * @param N
	 */
	virtual void allocate_particles(BaseParticle<number> **particles, int N) = 0;

	/**
	 * @brief Returns the maximum cut-off associated to the interaction
	 *
	 * @return Interaction cutoff
	 */
	virtual number get_rcut() { return _rcut; }

	/**
	 * @brief Check whether the initial configuration makes sense.
	 *
	 * Since this operation is interaction-dependent, each interaction must provide its own implementation.
	 * @param particles
	 * @param N
	 */
	virtual void check_input_sanity(BaseParticle<number> **particles, int N) = 0;

	/**
	 * @brief Computes the total interaction between particles p and q.
	 *
	 * If r is not given or NULL, it is computed from scratch. It can optionally update forces and torques exerted on p and q.
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return pair-interaction energy
	 */
	virtual number pair_interaction (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) = 0;
	
	/**
	 * @brief Computes the bonded part interaction between particles p and q.
	 *
	 * See {\@link pair_interaction} for a description of the parameters.
	 */
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) = 0;
	
	/**
	 * @brief Computed the non bonded part of the interaction between particles p and q.
	 *
	 * See {\@link pair_interaction} for a description of the parameters.
	 */
	virtual number pair_interaction_nonbonded (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) = 0;

	/**
	 * @brief Computes the requested term of the interaction energy between p and q.
	 *
	 * @param name identifier of the interaction method
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) = 0;

	/**
	 * @brief Returns the total potential energy of the system
	 *
	 * @param particles
	 * @param N
	 * @return
	 */
	virtual number get_system_energy(BaseParticle<number> **particles, int N) { return get_system_energy(particles, N, _box_side); }

	/**
	 * @brief Returns the total potential energy of the system, given a box size
	 *
	 * @param particles
	 * @param N
	 * @param box_side
	 * @return
	 */
	virtual number get_system_energy(BaseParticle<number> **particles, int N, number box_side);

	/**
	 * @brief return the total system energy according to a single term
	 * in the interaction potential
	 *
	 * @param name identifier of the interaction method to call for the
	 * pair interaction
	 * @param particles
	 * @param N
	 * @return the energy in question
	 */
	virtual number get_system_energy_term(int name, BaseParticle<number> **particles, int N) { return get_system_energy_term(name, particles, N, _box_side); }

	/**
	 * @brief Like get_system_energy_term, just that the box size can be specified in this case
	 *
	 * @param name
	 * @param particles
	 * @param N
	 * @param box_side
	 * @return the required energy contribution
	 */
	virtual number get_system_energy_term(int name, BaseParticle<number> **particles, int N, number box_side);

	/**
	 * @brief Read the topology of the interactions. The defaut
	 * method is empty.
	 *
	 * This function is set up so that the topology of the interactions
	 * are read in the code. The different interactions might require
	 * different formats, so it is implemented here. The default
	 *
	 * @param N number of particles as read previously. Used for
	 * checking consistency.
	 * @param N_strands Pointer to an int which will be filled with the
	 * number of "strands" (unbreakable clusters of connected particles)
	 * @param particles array of particles.
	 */
	virtual void read_topology (int N, int *N_strands, BaseParticle<number> **particles);

	/**
	 * @brief Returns the number of particles, as written in the topology.
	 *
	 * @return number of particles
	 */
	virtual int get_N_from_topology();

	/**
	 * @brief Returns a map which stores all the energy contributions due to the different
	 * energy terms.
	 *
	 * @return map of the energy contributions
	 */
	virtual map<int, number> get_system_energy_split(BaseParticle<number> **particles, int N) = 0;

	/**
	 * @brief Like get_system_energy_split, just that the box size can be
	 * specified in this case.
	 *
	 * @return map of the energy contributions
	 */
	virtual map<int, number> get_system_energy_split(BaseParticle<number> **particles, int N, number box_side) = 0;

	/**
	 * @brief Returns a list of potentially interaction neighbouts of particle p. This is VERY expensive, since
	 * it creates the cells and destroys them each time
	 *
	 * @param p particle
	 * @param particles
	 * @param N
	 * @param box_side
	 * @return an array storing a list of particles.
	 */
	virtual std::vector<BaseParticle<number>*> get_neighbours (BaseParticle<number> * p, BaseParticle<number> **particles, int N, number box_side) = 0;

	/**
	 * @brief Returns a list of potentially interacting pairs
	 *
	 * @param particles
	 * @param N
	 * @param box_side
	 * @return an array pairs of particle pointers.
	 */
	virtual std::vector< std::pair <BaseParticle<number> *, BaseParticle<number>*> > get_potential_interactions (BaseParticle<number> **particles, int N, number box_side) = 0;
	
	/**
	 * @brief Returns the state of the interaction
	 */
	bool get_is_infinite () { return _is_infinite; }
	
	/**
	 * @brief Sets _is_infinite
	 *
	 * @param arg bool 
	 */
	void set_is_infinite (bool arg) { _is_infinite = arg; }
	
	/**
	 * @brief overlap criterion. Used only in the generation of configurations. Can be overloaded.
	 *
	 * The default implementation does not allow two particles to have an
	 * interaction strength greater than 100. More subtle implementations
	 * may be needed in special cases.
	 *
	 * @param p
	 * @param q
	 * @param box_side
	 * @return bool whether the two particles overlap
	 */
	virtual bool generate_random_configuration_overlap (BaseParticle<number> * p, BaseParticle<number> *q, number box_side);
	
	
	/**
	 * @brief generation of initial configurations. Can be overloaded.
	 *
	 * The default function creates a random configuration in the most simple way possible:
	 * puts particles in one at a time, and using {\@link generate_random_configuration_overlap} 
	 * to test if two particles overlap.
	 *
	 * @param particles
	 * @param N
	 * @param box_side
	 */
	virtual void generate_random_configuration (BaseParticle<number> **particles, int N, number box_side);
};

template<typename number>
IBaseInteraction<number>::IBaseInteraction() {
	_energy_threshold = (number) 100.f;
}

template<typename number>
IBaseInteraction<number>::~IBaseInteraction() {
	_delete_cells();
}

template<typename number>
void IBaseInteraction<number>::get_settings(input_file &inp) {
	getInputString(&inp, "topology", this->_topology_filename, 1);
	getInputNumber (&inp, "energy_threshold", &_energy_threshold, 0);  
}

template<typename number>
int IBaseInteraction<number>::get_N_from_topology() {
	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	topology.getline(line, 512);
	topology.close();
	int ret;
	sscanf(line, "%d %*d\n", &ret);
	return ret;
}

template<typename number>
void IBaseInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;
	allocate_particles(particles, N);
	for (int i = 0; i < N; i ++) {
	   particles[i]->index = i;
	   particles[i]->type = 0;
	   particles[i]->strand_id = i;
	}
}

template<typename number>
void IBaseInteraction<number>::_create_cells(BaseParticle<number> **particles, int N, number box_side, bool init_all_neighs) {
	int N_cells = _cells_N_side*_cells_N_side*_cells_N_side;
	if(box_side != _last_box_side) {
		_delete_cells();

		_cells_N_side = (int) floor(box_side / _rcut + (number)0.1f);

		if(_cells_N_side < 3) _cells_N_side = 3;
		if(_cells_N_side > 300) _cells_N_side = 300;

		N_cells = _cells_N_side*_cells_N_side*_cells_N_side;

		_cells_head = new int[N_cells];
		// the allocation of the memory for the 27 neighbors of each cell is done later on only if
		// the cell is not empty (i.e. if _cells_head[i] != P_INVALID)
		_cells_neigh = new int *[N_cells];
		_cells_next = new int[N];
		_cells_index = new int[N];

		for(int i = 0; i < N_cells; i++) _cells_head[i] = P_INVALID;
	}
	// we set to P_INVALID only the cells which contained at least one particle in the previous _create_cells() call
	else for(std::set<int>::iterator it = _not_empty_cells.begin(); it != _not_empty_cells.end(); it++) _cells_head[*it] = P_INVALID;

	for(int i = 0; i < N; i++) _cells_next[i] = P_INVALID;

	_not_empty_cells.clear();
	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = particles[i];
		int cell_index = (int) ((p->pos.x / box_side - floor(p->pos.x / box_side)) * (1.f - FLT_EPSILON) * _cells_N_side);
		cell_index += _cells_N_side * ((int) ((p->pos.y / box_side - floor(p->pos.y / box_side)) * (1.f - FLT_EPSILON) * _cells_N_side));
		cell_index += _cells_N_side * _cells_N_side * ((int) ((p->pos.z / box_side - floor(p->pos.z / box_side)) * (1.f - FLT_EPSILON) * _cells_N_side));

		_not_empty_cells.insert(cell_index);

		int old_head = _cells_head[cell_index];
		_cells_head[cell_index] = i;
		_cells_index[i] = cell_index;
		assert(cell_index < N_cells);
		_cells_next[i] = old_head;
	}

	// if we have to initialize all the neighbours, we add all the cells to the set
	if(init_all_neighs) for(int i = 0; i < N_cells; i++) _not_empty_cells.insert(i);

	// we initialize only neighbours of the cells containing at least one particle
	for(std::set<int>::iterator it = _not_empty_cells.begin(); it != _not_empty_cells.end(); it++) {
		int i = *it;
		_cells_neigh[i] = new int[27];
		int ind[3], loop_ind[3], nneigh = 0;
		ind[0] = i % _cells_N_side;
		ind[1] = (i / _cells_N_side) % _cells_N_side;
		ind[2] = i / (_cells_N_side * _cells_N_side);
		for(int j = -1; j < 2; j++) {
			loop_ind[0] = (ind[0] + j + _cells_N_side) % _cells_N_side;
			for(int k = -1; k < 2; k++) {
				loop_ind[1] = (ind[1] + k + _cells_N_side) % _cells_N_side;
				for(int l = -1; l < 2; l++) {
					loop_ind[2] = (ind[2] + l + _cells_N_side) % _cells_N_side;
					int loop_index = (loop_ind[2] * _cells_N_side + loop_ind[1]) * _cells_N_side + loop_ind[0];
					_cells_neigh[i][nneigh] = loop_index;
					nneigh++;
				}
			}
		}
	}

	_last_box_side = box_side;
}

template<typename number>
void IBaseInteraction<number>::_delete_cell_neighs() {
	for(std::set<int>::iterator it = _not_empty_cells.begin(); it != _not_empty_cells.end(); it++) {
		delete[] _cells_neigh[*it];
	}
}

template<typename number>
void IBaseInteraction<number>::_delete_cells() {
	if(_cells_head != NULL)  {
		delete[] _cells_head;
		delete[] _cells_neigh;
		delete[] _cells_next;
		delete[] _cells_index;
	}
}

template<typename number>
number IBaseInteraction<number>::get_system_energy(BaseParticle<number> **particles, int N, number box_side) {
	_create_cells (particles, N, box_side);

	number energy = (number) 0.f;

	for (int i = 0; i < N; i ++) {
		BaseParticle<number> *p = particles[i];
		energy += pair_interaction_bonded(p, P_VIRTUAL);
		for(int c = 0; c < 27; c ++) {
			int j = _cells_head[_cells_neigh[_cells_index[p->index]][c]];
			while (j != P_INVALID) {
				BaseParticle<number> *q = particles[j];
				if(p->index > q->index) energy += pair_interaction_nonbonded(p, q);
				if(this->get_is_infinite()) {
					_delete_cell_neighs();
					return energy;
				}
				j = _cells_next[q->index];
			}
		}
	}

	_delete_cell_neighs();
	
	return energy;
}

template<typename number>
number IBaseInteraction<number>::get_system_energy_term(int name, BaseParticle<number> **particles, int N, number box_side) {
	_create_cells(particles, N, box_side);

	number energy = (number) 0.f;

	for (int i = 0; i < N; i ++) {
		BaseParticle<number> *p = particles[i];
		for(int c = 0; c < 27; c ++) {
			int j = _cells_head[_cells_neigh[_cells_index[p->index]][c]];
			while (j != P_INVALID) {
				BaseParticle<number> *q = particles[j];
				if(p->index > q->index) energy += pair_interaction_term(name, p, q);
				if(this->get_is_infinite()) {
					_delete_cell_neighs();
					return energy;
				}
				j = _cells_next[q->index];
			}
		}
	}

	_delete_cell_neighs();

	return energy;
}

template<typename number>
bool IBaseInteraction<number>::generate_random_configuration_overlap(BaseParticle<number> * p, BaseParticle<number> * q, number box_side) {
	LR_vector<number> dr = p->pos.minimum_image(q->pos, box_side);
	
	if (dr.norm() >= this->_sqr_rcut) return false;
	
	number energy = pair_interaction(p, q, &dr, false);

	// in case we create an overlap, we reset the interaction state
	this->set_is_infinite(false);
	
	// if energy too large, reject
	if (energy > _energy_threshold) return true;
	else return false;
}

template<typename number>
void IBaseInteraction<number>::generate_random_configuration(BaseParticle<number> **particles, int N, number box_side) {
	this->_create_cells(particles, N, box_side, true);
	
	
	for (int i = 0; i < _cells_N_side*_cells_N_side*_cells_N_side; i ++) _cells_head[i] = P_INVALID;
	for (int i = 0; i < N; i ++) _cells_next[i] = P_INVALID;
	for (int i = 0; i < N; i ++) _cells_index[i] = -1;
	
	LR_matrix<number> R;

	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = particles[i];
		
		bool inserted = false;
		int cell_index;
		do {
			p->pos = LR_vector<number> (drand48()*box_side, drand48()*box_side, drand48()*box_side);
			// random orientation
			p->orientation = Utils::get_random_rotation_matrix_from_angle<number> (M_PI);
			p->orientation.orthonormalize();
			p->orientationT = p->orientation.get_transpose();
			
			p->set_positions ();

			p->set_ext_potential (0, box_side);
		
			cell_index = (int) ((p->pos.x / box_side - floor(p->pos.x / box_side)) * (1.f - FLT_EPSILON) * this->_cells_N_side);
			cell_index += this->_cells_N_side * ((int) ((p->pos.y / box_side - floor(p->pos.y / box_side)) * (1.f - FLT_EPSILON) * this->_cells_N_side));
			cell_index += this->_cells_N_side * this->_cells_N_side * ((int) ((p->pos.z / box_side - floor(p->pos.z / box_side)) * (1.f - FLT_EPSILON) * this->_cells_N_side));

			inserted = true;
			for(int c = 0; c < 27 && inserted; c ++) {
				int j = this->_cells_head[this->_cells_neigh[cell_index][c]];
				while (j != P_INVALID && inserted) {
					BaseParticle<number> *q = particles[j];
					if (p != q && generate_random_configuration_overlap (p, q, box_side)) inserted = false;
					j = this->_cells_next[q->index];
				}
			}
			
			// we take into account the external potential
			if (p->ext_potential > _energy_threshold) {
				//if (inserted) printf ("rejecting because of ext. potential\n");
				inserted = false;
			}

		} while(!inserted);

		int old_head = this->_cells_head[cell_index];
		this->_cells_head[cell_index] = i;
		this->_cells_index[i] = cell_index;
		this->_cells_next[i] = old_head;

	}
	
	this->_delete_cell_neighs();
}


/**
 * @brief Base class for managing particle-particle interactions. It is an abstract class.
 *
 * This is the abstract class from which the other interaction classes inherit.
 * The children classes must implement all the public pair_interaction_*
 * methods.  We make use of the Curiously Recurring Template Pattern (see
 * http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) to make
 * it possible to call a method of any derived class from within this class
 * without knowing any detail of that derived class. Hence, we use the {\@link
 * _int_map} map to store function pointers of type child::*.
 * This class lacks a .cpp because we do not want to have template
 * instantiation of the children classes from within here.
 */
template <typename number, typename child>
class BaseInteraction : public IBaseInteraction<number> {
	typedef std::map<int, number (child::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)> interaction_map;

private:

protected:
	/**
	 * Map of function pointers of individual terms of the particle-particle interaction. Note that
	 * these pointers are of type child::* (i.e. methods defined within the child class).
	 */
	interaction_map _int_map;

	/**
	 * @brief Calls the right interaction-computing method, chosen according to its name.
	 *
	 * @param that instance of the child class that actually calls the right method
	 * @param name identifier of the interaction method
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return pair-interaction energy
	 */
	virtual number _pair_interaction_term_wrapper(child *that, int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

	/**
	 * @brief Build a mesh by using a function and its derivative.
	 *
	 * Only derived classes can call this function.
	 * @param that pointer to the class owning f and der
	 * @param f function
	 * @param der derivative of f
	 * @param pars pointer to a structure which contains function parameters
	 * @param npoints size of the mesh
	 * @param xlow the mesh is defined on a finite interval. This is the lower end
	 * @param xupp Upper end of the mesh interval
	 * @param m mesh to be built
	 */
	virtual void _build_mesh(child *that, number (child::*f)(number, void*), number (child::*der) (number, void*), void *pars, int npoints, number xlow, number xupp, Mesh<number> &m);
	inline number _query_mesh(number x, Mesh<number> &m);
	inline number _query_meshD(number x, Mesh<number> &m);

public:
	/**
	 * @brief Basic constructor. By default, it does not need anything.
	 */
	BaseInteraction();
	virtual ~BaseInteraction();


	/**
	 * @brief Like get_system_energy_split, just that the box size can be specified in this case
	 *
	 * @param particles
	 * @param N
	 * @return a map storing all the potential energy contributions.
	 */
	virtual map<int, number> get_system_energy_split(BaseParticle<number> **particles, int N) { return get_system_energy_split(particles, N, this->_box_side); }

	/**
	 * @brief Computes all the contributions to the total potential energy and returns them in a map.
	 *
	 * @param particles
	 * @param N
	 * @param box_side
	 * @return a map storing all the potential energy contributions.
	 */
	virtual map<int, number> get_system_energy_split(BaseParticle<number> **particles, int N, number box_side);

	virtual std::vector<BaseParticle<number>*> get_neighbours (BaseParticle<number> * p, BaseParticle<number> **particles, int N, number box_side);

	virtual std::vector<std::pair<BaseParticle<number> *, BaseParticle<number> *> > get_potential_interactions (BaseParticle<number> **particles, int N, number box_side);
	
};


template<typename number, typename child>
BaseInteraction<number, child>::BaseInteraction() : IBaseInteraction<number>() {
	this->_rcut = 2.5;
	this->_sqr_rcut = SQR(this->_rcut);
	this->_last_box_side = (number) 0;
	this->_cells_head = NULL;
	this->_cells_neigh = NULL;
	this->_cells_next = NULL;
	this->_cells_index = NULL;
	this->_is_infinite = false;
}

template<typename number, typename child>
BaseInteraction<number, child>::~BaseInteraction() {

}

template<typename number, typename child>
number BaseInteraction<number, child>::_pair_interaction_term_wrapper(child *that, int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r;
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
		r = &computed_r;
	}

	// _int_map.find(name) returns an iterator which points to the
	// requested function pointer. We dereference it, apply it to the
	// that, which is a pointer to an instance of a class which inherits
	// from BaseInteraction. We then pass in the required parameters and
	// return the interaction energy.
	typename interaction_map::iterator interaction = _int_map.find(name);
	if(interaction == _int_map.end()) {
		throw oxDNAException("%s, line %d: Interaction term '%d' not found", __FILE__, __LINE__, name);
	}

	return (that->*(interaction->second))(p, q, r, update_forces);
}

template<typename number, typename child>
inline number BaseInteraction<number, child>::_query_meshD(number x, Mesh<number> &m) {
	if (x < m.xlow) return m.B[0];
	if (x >= m.xupp) x = m.xupp - FLT_EPSILON;
	int i = (int) ((x - m.xlow)/m.delta);
	number dx = x - m.xlow - m.delta*i;
	return (m.B[i] + (2*dx*m.C[i] + 3*dx*dx*m.D[i]) * m.inv_sqr_delta);
}

template<typename number, typename child>
inline number BaseInteraction<number, child>::_query_mesh(number x, Mesh<number> &m) {
	if (x <= m.xlow) return m.A[0];
	if (x >= m.xupp) x = m.xupp - FLT_EPSILON;
	int i = (int) ((x - m.xlow)/m.delta);
	number dx = x - m.xlow - m.delta*i;
	return (m.A[i] + dx*(m.B[i] + dx*(m.C[i] + dx*m.D[i]) * m.inv_sqr_delta));
}

template<typename number, typename child>
void BaseInteraction<number, child>::_build_mesh(child *that, number (child::*f)(number, void*), number (child::*der) (number, void*), void * args, int npoints, number xlow, number xupp, Mesh<number> &m) {
    assert (xlow < xupp);
	int i;
	number x;

	m.init(npoints);

	number dx = (xupp - xlow) / (number) npoints;
	m.delta = dx;
	m.inv_sqr_delta = 1 / SQR(dx);
	m.xlow = xlow;
	m.xupp = xupp;

	number fx0, fx1, derx0, derx1;

	for(i = 0; i < npoints + 1; i++) {
		x = xlow + i * dx;

		fx0 = (that->*f)(x, args);
		fx1 = (that->*f)(x + dx, args);
		derx0 = (that->*der)(x, args);
		derx1 = (that->*der)(x + dx, args);

		m.A[i] = fx0;
		m.B[i] = derx0;
		m.D[i] = (2*(fx0 - fx1) + (derx0 + derx1)*dx) / dx;
		m.C[i] = (fx1 - fx0 + (-derx0 - m.D[i]) * dx);
    }
}

template<typename number, typename child>
map<int, number> BaseInteraction<number, child>::get_system_energy_split(BaseParticle<number> **particles, int N, number box_side) {
	this->_create_cells(particles, N, box_side);

	std::map<int, number> energy_map;

	for(typename interaction_map::iterator it = _int_map.begin(); it != _int_map.end(); it++) {
		int name = it->first;
		energy_map[name] = (number) 0.f;
	}

	for (int i = 0; i < N; i ++) {
		BaseParticle<number> *p = particles[i];
		for(int c = 0; c < 27; c ++) {
			int j = this->_cells_head[this->_cells_neigh[this->_cells_index[p->index]][c]];
			while (j != P_INVALID) {
				BaseParticle<number> *q = particles[j];
				if(q->index > p->index) {
					for(typename interaction_map::iterator it = _int_map.begin(); it != _int_map.end(); it++) {
						int name = it->first;
						energy_map[name] += this->pair_interaction_term(name, p, q);
					}
				}

				j = this->_cells_next[q->index];
			}
		}
	}

	this->_delete_cell_neighs();

	return energy_map;
}

template<typename number, typename child>
std::vector<BaseParticle<number> *> BaseInteraction<number, child>::get_neighbours(BaseParticle<number> *p, BaseParticle<number> **particles, int N, number box_side) {
	this->_create_cells(particles, N, box_side);

	std::vector <BaseParticle<number> *> ret;

	for(int c = 0; c < 27; c ++) {
		int j = this->_cells_head[this->_cells_neigh[this->_cells_index[p->index]][c]];
		while (j != P_INVALID) {
			BaseParticle<number> *q = particles[j];
			if (p != q) ret.push_back(q);
			j = this->_cells_next[q->index];
		}
	}

	this->_delete_cell_neighs();

	return ret;
}

template<typename number, typename child>
std::vector<std::pair<BaseParticle<number> *, BaseParticle<number> * > > BaseInteraction<number, child>::get_potential_interactions (BaseParticle<number> **particles, int N, number box_side) {
	this->_create_cells(particles, N, box_side);

	std::vector <std::pair<BaseParticle<number> *, BaseParticle<number> *> > ret;
	LR_vector<number> dr;

	for (int i = 0; i < N; i ++) {
		BaseParticle<number> *p = particles[i];
		for(int c = 0; c < 27; c ++) {
			int j = this->_cells_head[this->_cells_neigh[this->_cells_index[p->index]][c]];
			while (j != P_INVALID) {
				BaseParticle<number> *q = particles[j];
				if (p->index < q-> index) {
					//dr = q->pos.minimum_image(p->pos, this->_box_side);
					dr = q->pos.minimum_image(p->pos, box_side);
					if (dr.norm() < this->_sqr_rcut) ret.push_back(std::pair<BaseParticle<number>*, BaseParticle<number>*> (p, q));
				}
				j = this->_cells_next[q->index];
			}
		}
	}
	
	this->_delete_cell_neighs();

	return ret;
}
#endif /* BASE_INTERACTION_H */

