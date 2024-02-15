#include "CUDAPHBInteraction.cuh"


CUDAPHBInteraction::CUDAPHBInteraction() {
    _edge_compatible = true;
}

CUDAPHBInteraction::~CUDAPHBInteraction() {
}

void CUDAPHBInteraction::get_settings(input_file &inp) {
    PHBInteraction::get_settings(inp);
}

void CUDAPHBInteraction::get_cuda_settings(input_file &inp) {
    PHBInteraction::get_settings(inp);
}

void CUDAPHBInteraction::cuda_init(int N){
    CUDABaseInteraction::cuda_init(N);
    PHBInteraction::init();
}

void CUDAPHBInteraction::compute_forces() {
    
}