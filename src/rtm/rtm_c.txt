#ifndef RTM_C_H
#define RTM_C_H

#include "rtm.h"
#include "propagation.h"

void RTM_ForwardPropagation(rtm_t* r, int t);
void RTM_GetSourceSnapshots(rtm_t* r, int t);

void RTM_Backward_Propagation(rtm_t* r, int t);
void RTM_Accumulate_CrossCorrelation(rtm_t* r, int t);

#endif /* end of include guard: RTM_C_H */
