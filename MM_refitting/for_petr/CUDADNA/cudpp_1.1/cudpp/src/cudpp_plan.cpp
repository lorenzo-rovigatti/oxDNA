// -------------------------------------------------------------
// cuDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision: 3572$
// $Date: 2007-11-19 13:58:06 +0000 (Mon, 19 Nov 2007) $
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt
// in the root directory of this source distribution.
// ------------------------------------------------------------- 

#include "cudpp.h"
#include "cudpp_plan_manager.h"
#include "cudpp_scan.h"
#include "cudpp_segscan.h"
#include "cudpp_compact.h"
#include "cudpp_spmvmult.h"
#include "cudpp_radixsort.h"

#include <assert.h>

CUDPPPlanManager* CUDPPPlanManager::m_instance = NULL;

CUDPPResult validateOptions(CUDPPConfiguration config, size_t /*numElements*/, size_t numRows, size_t /*rowPitch*/)
{
    CUDPPResult ret = CUDPP_SUCCESS;
    if ((config.options & CUDPP_OPTION_BACKWARD) && (config.options & CUDPP_OPTION_FORWARD))
        ret = CUDPP_ERROR_ILLEGAL_CONFIGURATION;
    if ((config.options & CUDPP_OPTION_EXCLUSIVE) && (config.options & CUDPP_OPTION_INCLUSIVE))
        ret = CUDPP_ERROR_ILLEGAL_CONFIGURATION;

    if (config.algorithm == CUDPP_COMPACT && numRows > 1)
        ret = CUDPP_ERROR_ILLEGAL_CONFIGURATION; //!< @todo: add support for multi-row cudppCompact

    return ret;
}

/** @addtogroup publicInterface
  * @{
  */

/** @name Plan Interface
 * @{
 */


/** @brief Create a CUDPP plan 
  * 
  * A plan is a data structure containing state and intermediate storage space
  * that CUDPP uses to execute algorithms on data.  A plan is created by 
  * passing to cudppPlan() a CUDPPConfiguration that specifies the algorithm,
  * operator, datatype, and options.  The size of the data must also be passed
  * to cudppPlan(), in the \a numElements, \a numRows, and \a rowPitch 
  * arguments.  These sizes are used to allocate internal storage space at the
  * time the plan is created.  The CUDPP planner may use the sizes, options,
  * and information about the present hardware to choose optimal settings.
  *
  * Note that \a numElements is the maximum size of the array to be processed
  * with this plan.  That means that a plan may be re-used to process (for 
  * example, to sort or scan) smaller arrays.  
  * 
  * @param[out] planHandle A pointer to an opaque handle to the internal plan
  * @param[in]  config The configuration struct specifying algorithm and options
  * @param[in]  numElements The maximum number of elements to be processed
  * @param[in]  numRows The number of rows (for 2D operations) to be processed
  * @param[in]  rowPitch The pitch of the rows of input data, in elements
  */
CUDPP_DLL
CUDPPResult cudppPlan(CUDPPHandle        *planHandle, 
                      CUDPPConfiguration config, 
                      size_t             numElements, 
                      size_t             numRows, 
                      size_t             rowPitch)
{
    CUDPPResult result = CUDPP_SUCCESS;

    CUDPPPlan *plan;

    result = validateOptions(config, numElements, numRows, rowPitch);
    if (result != CUDPP_SUCCESS)
    {
        *planHandle = CUDPP_INVALID_HANDLE;
        return result;
    }

    switch (config.algorithm)
    {
    case CUDPP_SCAN:
        {
            plan = new CUDPPScanPlan(config, numElements, numRows, rowPitch);
            break;
        }
    case CUDPP_SORT_RADIX:
    //case CUDPP_SORT_RADIX_GLOBAL:
        {
            plan = new CUDPPRadixSortPlan(config, numElements);
            break;
        }
    //new rand plan
    case CUDPP_RAND_MD5:
        {
            plan = new CUDPPRandPlan(config, numElements);
            break;
        }
    case CUDPP_REDUCE:
    default:
        //! @todo: implement cudppReduce()
        return CUDPP_ERROR_ILLEGAL_CONFIGURATION; 
        break;
    }

    *planHandle = CUDPPPlanManager::AddPlan(plan);
    if (CUDPP_INVALID_HANDLE == *planHandle)
        return CUDPP_ERROR_UNKNOWN;
    else
        return CUDPP_SUCCESS;
}

/** @brief Destroy a CUDPP Plan
  *
  * Deletes the plan referred to by \a planHandle and all associated internal
  * storage.
  * 
  * @param[in] planHandle The CUDPPHandle to the plan to be destroyed
  */
CUDPP_DLL
CUDPPResult cudppDestroyPlan(CUDPPHandle planHandle)
{
    if (CUDPPPlanManager::RemovePlan(planHandle) == false)
        return CUDPP_ERROR_INVALID_HANDLE;
    else
        return CUDPP_SUCCESS;
}

/** @brief Plan base class constructor
  * 
  * @param[in]  config The configuration struct specifying algorithm and options
  * @param[in]  numElements The maximum number of elements to be processed
  * @param[in]  numRows The number of rows (for 2D operations) to be processed
  * @param[in]  rowPitch The pitch of the rows of input data, in elements
  */
CUDPPPlan::CUDPPPlan(CUDPPConfiguration config, 
                     size_t numElements, 
                     size_t numRows, 
                     size_t rowPitch)
: m_config(config),
  m_numElements(numElements),
  m_numRows(numRows),
  m_rowPitch(rowPitch)
{
}

/** @brief Scan Plan constructor
* 
* @param[in]  config The configuration struct specifying algorithm and options
* @param[in]  numElements The maximum number of elements to be scanned
* @param[in]  numRows The maximum number of rows (for 2D operations) to be scanned
* @param[in]  rowPitch The pitch of the rows of input data, in elements
*/
CUDPPScanPlan::CUDPPScanPlan(CUDPPConfiguration config, 
                             size_t numElements, 
                             size_t numRows, 
                             size_t rowPitch)
: CUDPPPlan(config, numElements, numRows, rowPitch),
  m_blockSums(0),
  m_rowPitches(0),
  m_numEltsAllocated(0),
  m_numRowsAllocated(0),
  m_numLevelsAllocated(0)
{
    allocScanStorage(this);
}

/** @brief CUDPP scan plan destructor */
CUDPPScanPlan::~CUDPPScanPlan()
{
    freeScanStorage(this);
}

/** @brief Sort Plan constructor
* 
* @param[in]  config The configuration struct specifying algorithm and options
* @param[in]  numElements The maximum number of elements to be sorted
*/
/*CUDPPSortPlan::CUDPPSortPlan(CUDPPConfiguration config, size_t numElements)
: CUDPPPlan(config, numElements, 1, 0),
  m_scanPlan(0),
  m_d_temp(0),
  m_d_tempAddress(0)
{
    CUDPPConfiguration scanConfig = 
    { 
      CUDPP_SCAN, 
      CUDPP_ADD, 
      CUDPP_UINT, 
      CUDPP_OPTION_FORWARD | CUDPP_OPTION_EXCLUSIVE 
    };

    //if (config.algorithm == CUDPP_SORT_RADIX_GLOBAL)
    {
        m_scanPlan = new CUDPPScanPlan(scanConfig, numElements, 1, 0);
    }

    allocSortStorage(this);
}*/

/** @brief Sort plan destructor */
/*CUDPPSortPlan::~CUDPPSortPlan()
{
    delete m_scanPlan;
    freeSortStorage(this);
}*/

CUDPPRadixSortPlan::CUDPPRadixSortPlan(CUDPPConfiguration config, size_t numElements)
: CUDPPPlan(config, numElements, 1, 0),
  m_scanPlan(0),
  m_tempKeys(0),    
  m_tempValues(0),
  m_counters(0),
  m_countersSum(0),
  m_blockOffsets(0) 
{
    size_t numBlocks2 = ((numElements % (SORT_CTA_SIZE * 2)) == 0) ?
            (numElements / (SORT_CTA_SIZE * 2)) : (numElements / (SORT_CTA_SIZE * 2) + 1);

    CUDPPConfiguration scanConfig = 
    { 
      CUDPP_SCAN, 
      CUDPP_ADD, 
      CUDPP_UINT, 
      CUDPP_OPTION_FORWARD | CUDPP_OPTION_EXCLUSIVE 
    };    

    if(m_config.options == CUDPP_OPTION_KEYS_ONLY)
        m_bKeysOnly = true;
    else
        m_bKeysOnly = false;

    m_scanPlan = new CUDPPScanPlan(scanConfig, numBlocks2*16, 1, 0);    
        
    allocRadixSortStorage(this); 
}

CUDPPRadixSortPlan::~CUDPPRadixSortPlan()
{
    delete m_scanPlan;
    freeRadixSortStorage(this);
}

/** @brief CUDPP Rand Plan Constructor
  * @param[in] config The configuration struct specifying options
  * @param[in] num_elements The number of elements to generate random bits for
  */
CUDPPRandPlan::CUDPPRandPlan(CUDPPConfiguration config, size_t num_elements) 
 : CUDPPPlan(config, num_elements, 1, 0),
   m_seed(0)
{
    
}

