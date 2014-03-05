#include "DNAInteraction_nomesh.h"

template<typename number>
DNAInteraction_nomesh<number>::DNAInteraction_nomesh() : DNAInteraction<number>() {

}

template<typename number>
DNAInteraction_nomesh<number>::~DNAInteraction_nomesh() {

}

template class DNAInteraction_nomesh<float>;
template class DNAInteraction_nomesh<double>;

