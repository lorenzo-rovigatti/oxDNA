// -------------------------------------------------------------
// cuDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision: 5632 $
// $Date: 2009-07-01 14:36:01 +1000 (Wed, 01 Jul 2009) $
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt in
// the root directory of this source distribution.
// ------------------------------------------------------------- 

/**
 * @file
 * cudpp.cpp
 *
 * @brief Main library source file.  Implements wrappers for public
 * interface.  
 * 
 * Main library source file.  Implements wrappers for public
 * interface.  These wrappers call application-level operators.
 * As this grows we may decide to partition into multiple source
 * files.
 */

/**
 * \defgroup publicInterface CUDPP Public Interface
 * The CUDA public interface comprises the functions, structs, and enums
 * defined in cudpp.h.  Public interface functions call functions in the
 * \link cudpp_app Application-Level\endlink interface. The public 
 * interface functions include Plan Interface functions and Algorithm
 * Interface functions.  Plan Inteface functions are used for creating
 * CUDPP Plan objects which contain configuration details, intermediate
 * storage space, and in the case of cudppSparseMatrix(), data.  The 
 * Algorithm Interface is the set of functions that do the real work 
 * of CUDPP, such as cudppScan() and cudppSparseMatrixVectorMultiply.
 *
 * @{
 */

/** @name Algorithm Interface
 * @{
 */

#include "cudpp.h"
#include "cudpp_plan_manager.h"
#include "cudpp_scan.h"
#include "cudpp_segscan.h"
#include "cudpp_compact.h"
#include "cudpp_spmvmult.h"
#include "cudpp_radixsort.h"
#include "cudpp_rand.h"

/**
 * @brief Performs a scan operation of numElements on its input in
 * GPU memory (d_in) and places the output in GPU memory
 * (d_out), with the scan parameters specified in the plan pointed to by
 * planHandle. 
 
 * The input to a scan operation is an input array, a binary associative 
 * operator (like + or max), and an identity element for that operator 
 * (+'s identity is 0). The output of scan is the same size as its input.
 * Informally, the output at each element is the result of operator
 * applied to each input that comes before it. For instance, the
 * output of sum-scan at each element is the sum of all the input
 * elements before that input.
 *
 * More formally, for associative operator
 * @htmlonly&oplus;@endhtmlonly@latexonly$\oplus$@endlatexonly,
 * <var>out<sub>i</sub></var> = <var>in<sub>0</sub></var>
 * @htmlonly&oplus;@endhtmlonly@latexonly$\oplus$@endlatexonly
 * <var>in<sub>1</sub></var>
 * @htmlonly&oplus;@endhtmlonly@latexonly$\oplus$@endlatexonly ...
 * @htmlonly&oplus;@endhtmlonly@latexonly$\oplus$@endlatexonly
 * <var>in<sub>i-1</sub></var>.
 * 
 * CUDPP supports "exclusive" and "inclusive" scans. For the ADD operator, 
 * an exclusive scan computes the sum of all input elements before the 
 * current element, while an inclusive scan computes the sum of all input 
 * elements up to and including the current element. 
 * 
 * Before calling scan, create an internal plan using cudppPlan().
 * 
 * After you are finished with the scan plan, clean up with cudppDestroyPlan(). 
 * 
 * @param[in] planHandle Handle to plan for this scan
 * @param[out] d_out output of scan, in GPU memory
 * @param[in] d_in input to scan, in GPU memory
 * @param[in] numElements number of elements to scan
 * 
 * @see cudppPlan, cudppDestroyPlan
 */
CUDPP_DLL
CUDPPResult cudppScan(CUDPPHandle planHandle,
                      void        *d_out, 
                      const void  *d_in, 
                      size_t      numElements)
{
    CUDPPScanPlan *plan = (CUDPPScanPlan*)CUDPPPlanManager::GetPlan(planHandle);
    if (plan != NULL)
    {
        cudppScanDispatch(d_out, d_in, numElements, 1, plan);
        return CUDPP_SUCCESS;
    }
    else
    {    
        return CUDPP_ERROR_UNKNOWN; //! @todo Return more specific errors
    }
}

/**
 * @brief Performs numRows parallel scan operations of numElements
 * each on its input (d_in) and places the output in d_out,
 * with the scan parameters set by config. Exactly like cudppScan 
 * except that it runs on multiple rows in parallel.
 * 
 * Note that to achieve good performance with cudppMultiScan one should
 * allocate the device arrays passed to it so that all rows are aligned
 * to the correct boundaries for the architecture the app is running on.
 * The easy way to do this is to use cudaMallocPitch() to allocate a 
 * 2D array on the device.  Use the \a rowPitch parameter to cudppPlan() 
 * to specify this pitch. The easiest way is to pass the device pitch 
 * returned by cudaMallocPitch to cudppPlan() via \a rowPitch.
 * 
 * @param[in] planHandle handle to CUDPPScanPlan
 * @param[out] d_out output of scan, in GPU memory
 * @param[in] d_in input to scan, in GPU memory
 * @param[in] numElements number of elements (per row) to scan
 * @param[in] numRows number of rows to scan in parallel
 * 
 * @see cudppScan, cudppPlan
 */
CUDPP_DLL
CUDPPResult cudppMultiScan(CUDPPHandle planHandle,
                            void       *d_out, 
                            const void *d_in, 
                            size_t     numElements,
                            size_t     numRows)
{
    CUDPPScanPlan *plan = (CUDPPScanPlan*)CUDPPPlanManager::GetPlan(planHandle);
    if (plan != NULL)
    {
        cudppScanDispatch(d_out, d_in, numElements, numRows, plan);
        return CUDPP_SUCCESS;
    }
    else
    {    
        return CUDPP_ERROR_UNKNOWN; //! @todo Return more specific errors
    }
}

/**
 * @brief Sorts key-value pairs or keys only
 * 
 * Takes as input an array of keys in GPU memory
 * (d_keys) and an optional array of corresponding values,
 * and outputs sorted arrays of keys and (optionally) values in place. 
 * Key-value and key-only sort is selected through the configuration of 
 * the plan, using the options CUDPP_OPTION_KEYS_ONLY and 
 * CUDPP_OPTION_KEY_VALUE_PAIRS.
 *
 * Supported key types are CUDPP_FLOAT and CUDPP_UINT.  Values can be
 * any 32-bit type (internally, values are treated only as a payload
 * and cast to unsigned int).
 *
 * @todo Determine if we need to provide an "out of place" sort interface.
 * 
 * @param[in] planHandle handle to CUDPPSortPlan
 * @param[out] d_keys keys by which key-value pairs will be sorted
 * @param[in] d_values values to be sorted
 * @param[in] keyBits the number of least significant bits in each element 
 *            of d_keys to sort by
 * @param[in] numElements number of elements in d_keys and d_values
 *
 * @see cudppPlan, CUDPPConfiguration, CUDPPAlgorithm
 */
CUDPP_DLL
CUDPPResult cudppSort(CUDPPHandle planHandle,
                      void        *d_keys,
                      void        *d_values,                      
                      int         keyBits,
                      size_t      numElements)
{
    CUDPPRadixSortPlan *plan = (CUDPPRadixSortPlan*)CUDPPPlanManager::GetPlan(planHandle);
    if (plan != NULL)
    {
        cudppRadixSortDispatch(d_keys, d_values, numElements, keyBits, plan);
        return CUDPP_SUCCESS;
    }
    else
    {
        return CUDPP_ERROR_UNKNOWN; //! @todo Return more specific errors.
    }
}

/**
 * @brief Rand puts \a numElements random 32-bit elements into \a d_out
 *
 
 * Outputs \a numElements random values to \a d_out. \a d_out must be of
 * type unsigned int, allocated in device memory.
 * 
 * The algorithm used for the random number generation is stored in \a planHandle.
 * Depending on the specification of the pseudo random number generator(PRNG),
 * the generator may have one or more seeds.  To set the seed, use cudppRandSeed().
 * 
 * @todo Currently only MD5 PRNG is supported.  We may provide more rand routines in 
 * the future.
 *
 * @param[in] planHandle Handle to plan for rand
 * @param[in] numElements number of elements in d_out.
 * @param[out] d_out output of rand, in GPU memory.  Should be an array of unsigned integers.
 *
 * @see cudppPlan, CUDPPConfiguration, CUDPPAlgorithm
 */
CUDPP_DLL
CUDPPResult cudppRand(CUDPPHandle planHandle,void * d_out, size_t numElements)
{
    CUDPPRandPlan * plan = (CUDPPRandPlan *) CUDPPPlanManager::GetPlan(planHandle);
    if(plan != NULL)
    {
        //dispatch the rand algorithm here
        cudppRandDispatch(d_out, numElements, plan);
        return CUDPP_SUCCESS;
    }
    else
        return CUDPP_ERROR_UNKNOWN; //! @todo Return more specific errors
}


/**@brief Sets the seed used for rand
 *
 * The seed is crucial to any random number generator as it allows a 
 * sequence of random numbers to be replicated.  Since there may be 
 * multiple different rand algorithms in CUDPP, cudppRandSeed 
 * uses \a planHandle to determine which seed to set.  Each rand 
 * algorithm has its own  unique set of seeds depending on what 
 * the algorithm needs.
 *
 * @param[in] planHandle the handle to the plan which specifies which rand seed to set
 * @param[in] seed the value which the internal cudpp seed will be set to
 */
CUDPP_DLL
CUDPPResult cudppRandSeed(const CUDPPHandle planHandle, unsigned int seed)
{
    CUDPPRandPlan * plan = (CUDPPRandPlan *) CUDPPPlanManager::GetPlan(planHandle);
    //switch on the plan to figure out which seed to update
    switch(plan->m_config.algorithm)
    {
    case CUDPP_RAND_MD5:
        plan->m_seed = seed;
        break;
    default:
        break;
    }

    return CUDPP_SUCCESS;
}//end cudppRandSeed

/** @} */ // end Algorithm Interface
/** @} */ // end of publicInterface group

// Leave this at the end of the file
// Local Variables:
// mode:c++
// c-file-style: "NVIDIA"
// End:

