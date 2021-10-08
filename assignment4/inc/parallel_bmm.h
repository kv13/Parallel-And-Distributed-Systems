#ifndef PARALLEL_BMM_H
#define PARALLEL_BMM_H


// ******* functions declaration ******* 
// ******* blocked bmm functions
void blocked_bmm_shared(coo_format *C, blocked_csr *blk_A, csr_format *AA, blocked_csc *blk_B, csc_format *BB);

// ******* filtered blocked bmm functions
void blocked_bmm_filtered_shared(coo_format *C, blocked_csr *blk_F,csr_format *FF,blocked_csr *blk_A, csr_format *AA, blocked_csc *blk_B, csc_format *BB);

#endif