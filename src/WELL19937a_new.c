/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/*                 A modified "maximally equidistributed" implementation         */
/*                 by Shin Harase, Hiroshima University.                         */
/* ***************************************************************************** */

#include <stdio.h>

#define W 32
#define R 624
#define DISCARD 31
#define MASKU (0xffffffffU>>(W-DISCARD))
#define MASKL (~MASKU)
#define M1 70
#define M2 179
#define M3 449

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define MAT1(v) v
#define MAT3POS(t,v) (v>>t)

#define V0            STATE[state_i]
#define VM1Over       STATE[state_i+M1-R]
#define VM1           STATE[state_i+M1]
#define VM2Over       STATE[state_i+M2-R]
#define VM2           STATE[state_i+M2]
#define VM3Over       STATE[state_i+M3-R]
#define VM3           STATE[state_i+M3]
#define VRm1          STATE[state_i-1]
#define VRm1Under     STATE[state_i+R-1]
#define VRm2          STATE[state_i-2]
#define VRm2Under     STATE[state_i+R-2]

#define newV0         STATE[state_i-1]
#define newV0Under    STATE[state_i-1+R]
#define newV1         STATE[state_i]
#define newVRm1       STATE[state_i-2]
#define newVRm1Under  STATE[state_i-2+R]

#define newVM2Over    STATE[state_i+M2-R+1]
#define newVM2        STATE[state_i+M2+1]

/*tempering paramater*/
#define BITMASK 0x41180000

/* #define FACT 2.32830643653869628906e-10 */
#define FACT 0x1.00000001P-32

static int state_i = 0;
static int func_num = 1;
static unsigned int STATE[R];
static unsigned int z0, z1, z2, y;
static double case_1 (void);
static double case_2 (void);
static double case_3 (void);
static double case_4 (void);
static double case_5 (void);
static double case_6 (void);
       double (*WELLRNG19937)(void);

void set_state (unsigned int *init_state, unsigned int *state_ind, unsigned int *func){
   int j;
   state_i = (int) *state_ind;

   switch(*func) {
      case 1:
        WELLRNG19937=case_1;
        func_num=1;
        break;
      case 2:
        WELLRNG19937=case_2;
        func_num=2;
        break;
      case 3:
        WELLRNG19937=case_3;
        func_num=3;
        break;
      case 4:
        WELLRNG19937=case_4;
        func_num=4;
        break;
      case 5:
        WELLRNG19937=case_5;
        func_num=5;
        break;
      case 6:
        WELLRNG19937=case_6;
        func_num=6;
        break;
      default:
        /* in case for inconsistent input */
        WELLRNG19937=case_1;
        func_num=1;
        break;
   }
   for (j = 0; j < R; j++)
     STATE[j] = init_state[j];
}

void get_state(unsigned int *init_state, unsigned int *state_ind, unsigned int *func)
{
   int j;
   *state_ind=state_i;
   *func=func_num; 
   for (j = 0; j<R; j++)
     init_state[j]=STATE[j];
}

static double case_1 (void){
   // state_i == 0
   z0 = (VRm1Under & MASKL) | (VRm2Under & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
   z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3);
   newV1      = z1 ^ z2;
   newV0Under = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i = R - 1;
   WELLRNG19937 = case_3;
   func_num=3;
   y = (STATE[state_i] ^ (newVM2Over & BITMASK));
   return ((double) y * FACT);
}

static double case_2 (void){
   // state_i == 1
   z0 = (VRm1 & MASKL) | (VRm2Under & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
   z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3);
   newV1 = z1 ^ z2;
   newV0 = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i = 0;
   WELLRNG19937 = case_1;
   func_num=1;
   y = (STATE[state_i] ^ (newVM2 & BITMASK));
   return ((double) y * FACT);
}

static double case_3 (void){
   // state_i+M1 >= R
   z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1Over);
   z2 = MAT3POS (9, VM2Over) ^ MAT0POS (1, VM3Over);
   newV1 = z1 ^ z2;
   newV0 = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i--;
   if (state_i + M1 < R) {
      WELLRNG19937 = case_5;
      func_num=5;
   }
   y = (STATE[state_i] ^ (newVM2Over & BITMASK));
   return ((double) y * FACT);
}

static double case_4 (void){
   // state_i+M3 >= R
   z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
   z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3Over);
   newV1 = z1 ^ z2;
   newV0 = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i--;
   if (state_i + M3 < R) {
      WELLRNG19937 = case_6;
      func_num=6;
   }
   y = (STATE[state_i] ^ (newVM2 & BITMASK));
   return ((double) y * FACT);
}

static double case_5 (void){
   // state_i+M2 >= R
   z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
   z2 = MAT3POS (9, VM2Over) ^ MAT0POS (1, VM3Over);
   newV1 = z1 ^ z2;
   newV0 = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i--;
   if (state_i + M2 < R) {
      WELLRNG19937 = case_4;
      func_num=4;
   }
   y = (STATE[state_i] ^ (newVM2Over & BITMASK));
   return ((double) y * FACT);
}

static double case_6 (void){
   // 2 <= state_i <= (R - M3 - 1)
   z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
   z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
   z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3);
   newV1 = z1 ^ z2;
   newV0 = MAT1 (z0) ^ MAT0NEG (-9, z1) ^ MAT0NEG (-21, z2) ^ MAT0POS (21, newV1);
   state_i--;
   if (state_i == 1) {
      WELLRNG19937 = case_2;
      func_num=2;
   }
   y = (STATE[state_i] ^ (newVM2 & BITMASK));
   return ((double) y * FACT);
}
