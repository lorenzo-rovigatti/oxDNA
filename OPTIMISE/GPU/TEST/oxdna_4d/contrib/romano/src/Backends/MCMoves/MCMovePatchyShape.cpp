/**
 * @file    MCMovePatchyShape.cpp
 * @date    30/apr/2014
 * @author  flavio
 *
 *
 */
#include "MCMovePatchyShape.h"

/// traslation
template<typename number>
MCMovePatchyShape<number>::MCMovePatchyShape (){
	pos_old = LR_vector<number> (0., 0., 0.);

	_verlet_skin = -1.f;
}

template<typename number>
MCMovePatchyShape<number>::~MCMovePatchyShape () {

}


template<typename number>
void MCMovePatchyShape<number>::init () {
	BaseMove<number>::init();
	PatchyShapeInteraction<number> *interaction = dynamic_cast<PatchyShapeInteraction<number> *  >( this->_Info->interaction );
	interaction->_init_patchy_locks();
	//interaction->check_patchy_locks();
}

template<typename number>
void MCMovePatchyShape<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta_translation", &_delta, 1);
	getInputNumber (&inp, "delta_rotation", &_delta_rotation, 1);

	OX_LOG(Logger::LOG_INFO, "(MCMovePatchyShape.cpp) MCMoveShapeParticle initiated with T %g, delta_translation %g, delta_rotation %g, prob: %g", this->_T, _delta,_delta_rotation, this->prob);
	
	std::string tmps;
	if (getInputString (&sim_inp, "list_type", tmps, 0) == KEY_FOUND) {
		if (!tmps.compare ("verlet")) {
			getInputNumber (&sim_inp, "verlet_skin", &_verlet_skin, 1);
		}
	}
}


struct lock
{
	int pid;
	int ppatchid;
	int qid;
	int qpatchid;

	lock(int _pid, int _qid, int plockid, int qlockid) : pid(_pid),  ppatchid(plockid), qid(_qid), qpatchid(qlockid) {}

	bool operator == (lock &b)
	{
		return bool( (pid == b.pid && ppatchid == b.ppatchid && qid == b.qid && qpatchid == b.qpatchid ) || (qid == b.pid && qpatchid == b.ppatchid && pid == b.qid && ppatchid == b.qpatchid ) );
	}
};

template<typename number>
void MCMovePatchyShape<number>::apply (llint curr_step) {

	// we increase the attempted count
	this->_attempted += 1;

	// we select the particle to translate
	int pi = (int) (drand48() * (*this->_Info->N));

	//cout << "GGB " << this->_Info->particles << endl;;

	PatchyShapeParticle<number> *p = static_cast< PatchyShapeParticle<number> *>(this->_Info->particles[pi]);

	//printf("Moving particle %d with pos %f \n",p->index,p->pos[0]);

	pos_old = p->pos; 
	_orientation_old = p->orientation;
	_orientationT_old = p->orientationT;

	// compute the energy before the move
	number delta_E = -this->particle_energy(p);
	p->set_ext_potential (curr_step, this->_Info->box);
	number delta_E_ext = -p->ext_potential;

	//printf("Before suggesting move: delta_E = %f, pos = %g %g %g \n",delta_E,p->pos.x,p->pos.y,p->pos.z);
	//printf("Using delta_trans: %g  delta_rot %g" ,_delta,_delta_rotation);
	// perform the move
	if(drand48() < 0.5) // translation
	{
	 p->pos.x += 2. * (drand48() - (number)0.5f) * _delta;
	 p->pos.y += 2. * (drand48() - (number)0.5f) * _delta;
	 p->pos.z += 2. * (drand48() - (number)0.5f) * _delta;
	}
	else { //rotation

		number t = drand48() * _delta_rotation;
		LR_vector<number> axis = Utils::get_random_vector<number>();

		number sintheta = sin(t);
		number costheta = cos(t);
		number olcos = ((number)1.) - costheta;

		number xyo = axis.x * axis.y * olcos;
		number xzo = axis.x * axis.z * olcos;
		number yzo = axis.y * axis.z * olcos;
		number xsin = axis.x * sintheta;
		number ysin = axis.y * sintheta;
		number zsin = axis.z * sintheta;

		LR_matrix<number> R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin,
				xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin,
				xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

		p->orientation = p->orientation * R;
		p->orientationT = p->orientation.get_transpose();
		p->set_positions();

	}
	// update lists
	this->_Info->lists->single_update(p);
	if(!this->_Info->lists->is_updated()) {
		this->_Info->lists->global_update();
	}

	//vector< Patch<number>  > broken_locks;
	vector<lock> locks_to_restore; //this is a list of locks that will be rejected
    vector<lock> locks_to_break;
	PatchyShapeInteraction<number> *interaction = static_cast<PatchyShapeInteraction<number> *  >(this->_Info->interaction);

	//printf("After suggesting mov and fucking: , delta_E = %f, pos = %g %g %g \n",delta_E,p->pos.x,p->pos.y,p->pos.z);
	for(int i = 0; i < p->N_patches; i++) {
		if(p->patches[i].is_locked()) {
			int qid=-1, qpatch=-1;
			p->patches[i].get_lock(qid,qpatch);
			assert(qid >= 0);  assert(qpatch >= 0 && qpatch < 6);
			PatchyShapeParticle<number> *q = static_cast< PatchyShapeParticle<number> *>(this->_Info->particles[qid]);
			LR_vector<number> r = this->_Info->box->min_image(p,q);
			number new_ene = interaction->just_two_patch_interaction(p,q,i,qpatch,&r);
			//printf("I am particle %d, patch %d, locked to %d %d (new energy %g), cutoff %f\n",p->index,i,qid,qpatch, new_ene,interaction->get_patch_cutoff_energy() );
			//if ( new_ene < -0.1 ) throw oxDNAException("it actually happens");
			if ( ! (new_ene < interaction->get_patch_cutoff_energy()) ) //we break the lock
			{

				//printf("I am particle %d breaking lock %d with %d (%d) \n",p->index,i,qid,qpatch);
				// broken_locks.push_back(p->patches[i]);
				p->patches[i].unlock();
				q->patches[qpatch].unlock();
				locks_to_restore.push_back(  lock(p->index,q->index,i,qpatch));
				//throw oxDNAException("here...");
			}
		}
	}

	// energy after the move
	delta_E += this->particle_energy(p);
	p->set_ext_potential(curr_step, this->_Info->box);
	delta_E_ext += p->ext_potential;

	number delta_E_newlocks = 0.f;

	for(vector<lock>::iterator i = locks_to_restore.begin(); i != locks_to_restore.end(); ++i)
	{
		//i->locked_to(qid,qpatch);
		// we use the locks_to_restore array just to get the indexes 
		// of the particles (and patches) to check
		int qid = i->qid;
		int qpatch = i->qpatchid;

		PatchyShapeParticle<number> *q = static_cast< PatchyShapeParticle<number> *>(this->_Info->particles[qid]);
		std::vector<BaseParticle<number> *> neighs = this->_Info->lists->get_neigh_list(q);
		for(unsigned int n = 0; n < neighs.size(); n++) {
			PatchyShapeParticle<number> *qq =   static_cast< PatchyShapeParticle<number> *>(neighs[n]);
			LR_vector<number> r = this->_Info->box->min_image(q,qq);
			for(int qqpatch = 0; qqpatch < qq->N_patches; qqpatch++)
			{
				number new_ene = interaction->just_two_patch_interaction(q,qq,qpatch,qqpatch,&r);
				if(new_ene < interaction->get_patch_cutoff_energy())
				{
					q->patches[qpatch].set_lock(qq->index,qqpatch);
					qq->patches[qqpatch].set_lock(qid,qpatch);
					//make sure the lock is not counted before
					lock newlock =  lock(q->index,qq->index,qpatch,qqpatch);
					bool known = false;
					for(vector<lock>::iterator it = locks_to_break.begin(); it != locks_to_break.end(); it++)
					{
						if(*it == newlock)
						{
							known = true;
							break;
						}
					}
					if (! known)
					{
					  delta_E_newlocks += new_ene;
					  locks_to_break.push_back(newlock);
					}
				}
			}
		}
	}

	//printf("After suggesting move: delta_E_newlocks: %f , delta_E = %f, pos = %g %g %g \n",delta_E_newlocks,delta_E,p->pos.x,p->pos.y,p->pos.z);
	// accept or reject?
	if (this->_Info->interaction->get_is_infinite() == false && ((delta_E + delta_E_ext + delta_E_newlocks) < 0 || exp(-(delta_E + delta_E_ext+delta_E_newlocks) / this->_T) > drand48() )) {
		// move accepted
        //printf("Move accepted\n");
		this->_accepted ++;
		//We need to make new locks for p
		std::vector<BaseParticle<number> *> neighs = this->_Info->lists->get_neigh_list(p);
		for(unsigned int n = 0; n < neighs.size(); n++) {
			PatchyShapeParticle<number> *qq = static_cast< PatchyShapeParticle<number> *>(neighs[n]);
			LR_vector<number> r = this->_Info->box->min_image(p,qq);
			for(int ppatch = 0; ppatch < p->N_patches; ppatch++ )
			{
				for(int qqpatch = 0; qqpatch < qq->N_patches; qqpatch++)
				{
					number new_ene = interaction->just_two_patch_interaction(p,qq,ppatch,qqpatch,&r);
					if(new_ene < interaction->get_patch_cutoff_energy())
					{
						p->patches[ppatch].set_lock(qq->index,qqpatch);
						qq->patches[qqpatch].set_lock(p->index,ppatch);
						// printf("setting lock %d (%d) %d (%d) \n",p->index,ppatch,qq->index,qqpatch);

					}
				}
			}
		}

		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			_delta *= this->_acc_fact;
			if (_verlet_skin > 0. && _delta > _verlet_skin * 0.8660254037844386 / 2.) {
				_delta = _verlet_skin *_verlet_skin * 0.8660254037844386 / 2.; //  0.8660254 .. ==sqrt(3.) / 2.
			}
		}
	}
	else {
		//printf("Move rejected\n");
		//rejected; We need to fix locks that were broken or created
		for (vector<lock>::iterator i = locks_to_break.begin(); i != locks_to_break.end(); ++i)
		{
			PatchyShapeParticle<number> *pp = static_cast< PatchyShapeParticle<number> *>(this->_Info->particles[i->pid]);
			PatchyShapeParticle<number> *qq = static_cast< PatchyShapeParticle<number> *>(this->_Info->particles[i->qid]);
			pp->patches[i->ppatchid].unlock();
			qq->patches[i->qpatchid].unlock();
			//printf("unsetting lock %d %d %d %d \n",pp->index,i->ppatchid,qq->index,i->qpatchid);
		}
		for (vector<lock>::iterator i = locks_to_restore.begin(); i != locks_to_restore.end(); ++i)
		{
			PatchyShapeParticle<number> *pp = static_cast< PatchyShapeParticle<number> *>(this->_Info->particles[i->pid]);
			PatchyShapeParticle<number> *qq = static_cast< PatchyShapeParticle<number> *>(this->_Info->particles[i->qid]);
			pp->patches[i->ppatchid].set_lock(i->qid,i->qpatchid);
			qq->patches[i->qpatchid].set_lock(i->pid,i->ppatchid);
			//printf("restoring lock %d %d %d %d \n",pp->index,i->ppatchid,qq->index,i->qpatchid);
		}

		this->_Info->particles[pi]->pos = pos_old;
		p->orientation = _orientation_old;
		p->orientationT = _orientationT_old;
		p->set_positions();

		this->_Info->lists->single_update(p);
		this->_Info->interaction->set_is_infinite(false);

		if (curr_step < this->_equilibration_steps && this->_adjust_moves) _delta /= this->_rej_fact;
	}

	//perform check
	if(3>4 && curr_step % 10000  == 0)
	{
		interaction->check_patchy_locks();
		//printf("all good\n");
	}
	return;
}

template class MCMovePatchyShape<float>;
template class MCMovePatchyShape<double>;
