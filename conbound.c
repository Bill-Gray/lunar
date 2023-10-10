#include <string.h>
#include <stdint.h>
#include "watdefs.h"
#include "afuncs.h"

/* See 'constbnd.c' in the 'constbnd' repository for a description of
how this all works.  Basically,  this is an highly optimized algorithm
for determining the constellation in which a given RA/dec falls.  The
algorithm (with some minor differences) is also described at

https://iopscience.iop.org/article/10.1086/132034/pdf

   Note that the constellation boundaries fall along RA/dec boundaries
in epoch B1875.  That's emphasized in the '*_degrees_1875' names of the
parameters passed to constell_from_ra_dec( ).

   Compile with

gcc -Wall -Wextra -pedantic -Werror -DTEST_MAIN -o conbound conbound.c

   for a small program to test/demonstrate this function.

   Three-letter abbreviations for 88 constellations = 264-byte string : */
const char *constell_names =
        "AndAntApsAqlAqrAraAriAurBooCaeCamCapCarCasCenCepCetChaCirCMaCMiCnc"
        "ColComCrACrBCrtCruCrvCVnCygDelDorDraEquEriForGemGruHerHorHyaHyiInd"
        "LacLeoLepLibLMiLupLynLyrMenMicMonMusNorOctOphOriPavPegPerPhePicPsA"
        "PscPupPyxRetSclScoSctSerSexSgeSgrTauTelTrATriTucUMaUMiVelVirVolVul";

static const struct
   {
   int32_t ra_spd;   /* ra1 in 17 bits, south polar distance in 14 bits */
   uint16_t d_ra;    /* ra2 - ra1 */
   char constell_idx;      /* from 0 to 87 */
   } bounds[324] = {
    { 0x53714370, 0x7e90, 15 },   /* Cep */
    { 0x52bc7080, 0x5b68, 10 },   /* Cam */
    { 0x52952750, 0x9ab0, 15 },   /* Cep */
    { 0x5280fd20, 0x2a30, 33 },   /* Dra */
    { 0x52084650, 0x8598, 10 },   /* Cam */
    { 0x50a080e8, 0x1518, 33 },   /* Dra */
    { 0x4fb11b98, 0x673e, 15 },   /* Cep */
    { 0x4fb0f618, 0x2580, 33 },   /* Dra */
    { 0x4fb0bf04, 0x3714, 83 },   /* UMi */
    { 0x4fb080e8, 0x20d0, 33 },   /* Dra */
    { 0x4fb03156, 0x4f92, 10 },   /* Cam */
    { 0x4e48b6d0, 0x3f48, 83 },   /* UMi */
    { 0x4e4880e8, 0x35e8, 33 },   /* Dra */
    { 0x4e48300c, 0x50dc, 10 },   /* Cam */
    { 0x4e4804b0, 0x2b5c, 13 },   /* Cas */
    { 0x4d58e880, 0x3a20, 33 },   /* Dra */
    { 0x4ca47008, 0x2f58, 82 },   /* UMa */
    { 0x4b00dc50, 0x4650, 33 },   /* Dra */
    { 0x4b009f60, 0x2580, 33 },   /* Dra */
    { 0x4a102b98, 0x4470, 10 },   /* Cam */
    { 0x49991f1c, 0x3714, 15 },   /* Cep */
    { 0x495c7008, 0x38b8, 82 },   /* UMa */
    { 0x49214ba4, 0x3174, 13 },   /* Cas */
    { 0x4920a8c0, 0x765c, 33 },   /* Dra */
    { 0x48307008, 0x4dd0, 82 },   /* UMa */
    { 0x47b945c8, 0x3750, 13 },   /* Cas */
    { 0x47b87008, 0x5ab4, 82 },   /* UMa */
    { 0x474055c8, 0x0ca8, 50 },   /* Lyn */
    { 0x47051940, 0x2c88, 15 },   /* Cep */
    { 0x46bf20cc, 0x00e4, 30 },   /* Cyg */
    { 0x465055c8, 0x2094, 50 },   /* Lyn */
    { 0x461515f8, 0x0bb8, 30 },   /* Cyg */
    { 0x45e34190, 0x3b88, 13 },   /* Cas */
    { 0x459c1ad6, 0x0762, 62 },   /* Per */
    { 0x4561110c, 0x10a4, 30 },   /* Cyg */
    { 0x452417e8, 0x0a50, 62 },   /* Per */
    { 0x44e817e8, 0x14a0, 62 },   /* Per */
    { 0x448f39d4, 0x07bc, 44 },   /* Lac */
    { 0x44704650, 0x0f78,  7 },   /* Aur */
    { 0x44350c5c, 0x1554, 30 },   /* Cyg */
    { 0x4434c558, 0x111c,  8 },   /* Boo */
    { 0x43f93740, 0x0a50, 44 },   /* Lac */
    { 0x43f817e8, 0x16f8, 62 },   /* Per */
    { 0x43e50c5c, 0x288c, 30 },   /* Cyg */
    { 0x43804650, 0x1518,  7 },   /* Aur */
    { 0x43801338, 0x1ba8, 62 },   /* Per */
    { 0x4308c558, 0x1824,  8 },   /* Boo */
    { 0x4308a9ec, 0x13ec, 29 },   /* CVn */
    { 0x42eb34e8, 0x0ca8, 44 },   /* Lac */
    { 0x42cd4190, 0x0690,  0 },   /* And */
    { 0x42cc41fa, 0x196e,  7 },   /* Aur */
    { 0x42cc1338, 0x2ec2, 62 },   /* Per */
    { 0x4254dd7c, 0x1194, 39 },   /* Her */
    { 0x41dcdd7c, 0x22ec, 39 },   /* Her */
    { 0x41dc1cb6, 0x06ae,  0 },   /* And */
    { 0x41a14190, 0x0a14,  0 },   /* And */
    { 0x41a041fa, 0x1da6,  7 },   /* Aur */
    { 0x41a00fb4, 0x07bc,  0 },   /* And */
    { 0x40eca9ec, 0x1b6c, 29 },   /* CVn */
    { 0x40b14190, 0x1248,  0 },   /* And */
    { 0x40b00c30, 0x0b40,  0 },   /* And */
    { 0x4074ff96, 0x0df2, 51 },   /* Lyr */
    { 0x40385fa0, 0x2148, 50 },   /* Lyn */
    { 0x40380c30, 0x1734,  0 },   /* And */
    { 0x3fc14190, 0x3354,  0 },   /* And */
    { 0x3f48a8c0, 0x1c98, 29 },   /* CVn */
    { 0x3f0c41fa, 0x259e,  7 },   /* Aur */
    { 0x3ed13416, 0x0d7a, 44 },   /* Lac */
    { 0x3eb3339e, 0x0df2, 44 },   /* Lac */
    { 0x3e94ff96, 0x113a, 51 },   /* Lyr */
    { 0x3de086c4, 0x0834, 48 },   /* LMi */
    { 0x3de06798, 0x1f2c, 50 },   /* Lyn */
    { 0x3cf0d908, 0x0ca8, 25 },   /* CrB */
    { 0x3cf086c4, 0x10e0, 48 },   /* LMi */
    { 0x3cd28214, 0x1590, 48 },   /* LMi */
    { 0x3b6a1c20, 0x07f8, 80 },   /* Tri */
    { 0x3b4d103a, 0x2364, 30 },   /* Cyg */
    { 0x3b1131a0, 0x03c0, 61 },   /* Peg */
    { 0x3b103f48, 0x2850,  7 },   /* Aur */
    { 0x3ad45be0, 0x111c, 37 },   /* Gem */
    { 0x3a9931a0, 0x0f3c, 61 },   /* Peg */
    { 0x3a9813ce, 0x104a, 80 },   /* Tri */
    { 0x3a5d31a0, 0x18d8, 61 },   /* Peg */
    { 0x3a20a8c0, 0x04b0, 23 },   /* Com */
    { 0x3a208214, 0x189c, 48 },   /* LMi */
    { 0x3a2013ce, 0x1266, 80 },   /* Tri */
    { 0x39e48214, 0x08e8, 45 },   /* Leo */
    { 0x39e47080, 0x1194, 21 },   /* Cnc */
    { 0x39e45be0, 0x14a0, 37 },   /* Gem */
    { 0x39a8d584, 0x102c, 25 },   /* CrB */
    { 0x39a80a14, 0x09ba, 66 },   /* Psc */
    { 0x393b31a0, 0x1c5c, 61 },   /* Peg */
    { 0x3930a8c0, 0x1194, 23 },   /* Com */
    { 0x38e131a0, 0x1fe0, 61 },   /* Peg */
    { 0x389ac44a, 0x113a,  8 },   /* Boo */
    { 0x38902f58, 0x0ff0, 77 },   /* Tau */
    { 0x389021fc, 0x0d5c,  6 },   /* Ari */
    { 0x38410ed2, 0x22ce, 30 },   /* Cyg */
    { 0x3840e5b0, 0x1c98, 39 },   /* Her */
    { 0x38402f58, 0x1374, 77 },   /* Tau */
    { 0x37c91490, 0x1194, 87 },   /* Vul */
    { 0x37c8a6e0, 0x1374, 23 },   /* Com */
    { 0x37c89ab0, 0x0c30, 45 },   /* Leo */
    { 0x378cbdd8, 0x17ac,  8 },   /* Boo */
    { 0x378ca6e0, 0x16f8, 23 },   /* Com */
    { 0x378c8214, 0x1194, 45 },   /* Leo */
    { 0x378c2f58, 0x2364, 77 },   /* Tau */
    { 0x37512d2c, 0x2544, 61 },   /* Peg */
    { 0x37511490, 0x189c, 87 },   /* Vul */
    { 0x37506edc, 0x1338, 21 },   /* Cnc */
    { 0x375052bc, 0x1c20, 37 },   /* Gem */
    { 0x37500a14, 0x0d5c, 66 },   /* Psc */
    { 0x37150ed2, 0x1e5a, 87 },   /* Vul */
    { 0x36f61af4, 0x1464,  6 },   /* Ari */
    { 0x36d8e358, 0x1ef0, 39 },   /* Her */
    { 0x3660e178, 0x27d8, 39 },   /* Her */
    { 0x3660d41c, 0x0d5c, 73 },   /* Ser */
    { 0x36250950, 0x23dc, 87 },   /* Vul */
    { 0x3624972c, 0x0fb4, 45 },   /* Leo */
    { 0x35e81770, 0x17e8,  6 },   /* Ari */
    { 0x355200f0, 0x0b04,  0 },   /* And */
    { 0x35352ad4, 0x279c, 61 },   /* Peg */
    { 0x35348214, 0x24cc, 45 },   /* Leo */
    { 0x34e45028, 0x0294, 59 },   /* Ori */
    { 0x34812ad4, 0x28aa, 61 },   /* Peg */
    { 0x3480dfd4, 0x297c, 39 },   /* Her */
    { 0x34445028, 0x0744, 59 },   /* Ori */
    { 0x342716e8, 0x05dc, 75 },   /* Sge */
    { 0x34130950, 0x0564, 75 },   /* Sge */
    { 0x340801fe, 0x1572, 66 },   /* Psc */
    { 0x33cd1cc4, 0x0474, 31 },   /* Del */
    { 0x33906dce, 0x1446, 21 },   /* Cnc */
    { 0x33552804, 0x2b7a, 61 },   /* Peg */
    { 0x33551cc4, 0x0b40, 31 },   /* Del */
    { 0x332d0950, 0x1374, 75 },   /* Sge */
    { 0x33182e2c, 0x21fc, 77 },   /* Tau */
    { 0x32dd0950, 0x01e0,  3 },   /* Aql */
    { 0x32a02e2c, 0x22ec, 77 },   /* Tau */
    { 0x32645118, 0x079e, 59 },   /* Ori */
    { 0x31c50950, 0x0d98,  3 },   /* Aql */
    { 0x31b0d41c, 0x0e10, 73 },   /* Ser */
    { 0x31b045d8, 0x0528, 59 },   /* Ori */
    { 0x31931b3e, 0x0cc6, 31 },   /* Del */
    { 0x31930950, 0x11ee,  3 },   /* Aql */
    { 0x317440ec, 0x0dd4, 59 },   /* Ori */
    { 0x3138b478, 0x0960, 85 },   /* Vir */
    { 0x30e8f294, 0x0e10, 58 },   /* Oph */
    { 0x30c0a6e0, 0x16f8, 85 },   /* Vir */
    { 0x30846978, 0x0456, 20 },   /* CMi */
    { 0x3034eb8c, 0x1518, 58 },   /* Oph */
    { 0x300d28f4, 0x030c, 34 },   /* Equ */
    { 0x300c6270, 0x0b5e, 20 },   /* CMi */
    { 0x300c40ec, 0x17ca, 59 },   /* Ori */
    { 0x300c0000, 0x1770, 66 },   /* Psc */
    { 0x2fd10670, 0x14ce,  3 },   /* Aql */
    { 0x2fd0eb8c, 0x1ae4, 58 },   /* Oph */
    { 0x2fd058b6, 0x08ca, 54 },   /* Mon */
    { 0x2fbd258e, 0x0672, 34 },   /* Equ */
    { 0x2f58a1f4, 0x1be4, 85 },   /* Vir */
    { 0x2ee14f28, 0x19c8, 66 },   /* Psc */
    { 0x2ee06270, 0x0d02, 20 },   /* CMi */
    { 0x2ee057c6, 0x0aaa, 54 },   /* Mon */
    { 0x2ed74f28, 0x1e78, 66 },   /* Psc */
    { 0x2ed61c20, 0x120c, 16 },   /* Cet */
    { 0x2e2d0670, 0x1708,  3 },   /* Aql */
    { 0x2df0a1f4, 0x3228, 85 },   /* Vir */
    { 0x2db53fec, 0x2db4, 66 },   /* Psc */
    { 0x2d7886c4, 0x1068, 74 },   /* Sex */
    { 0x2d7871ac, 0x1518, 41 },   /* Hya */
    { 0x2d786270, 0x0f3c, 20 },   /* CMi */
    { 0x2d1f00a4, 0x08ac, 73 },   /* Ser */
    { 0x2d0124f8, 0x0708, 34 },   /* Equ */
    { 0x2cc457c6, 0x0ae6, 54 },   /* Mon */
    { 0x2c4ceb8c, 0x178e, 58 },   /* Oph */
    { 0x2c10e4c0, 0x1e5a, 58 },   /* Oph */
    { 0x2c10d41c, 0x10a4, 73 },   /* Ser */
    { 0x2b9900a4, 0x08ac, 73 },   /* Ser */
    { 0x2b7b2de0, 0x02d0,  4 },   /* Aqr */
    { 0x2b213560, 0x0a8c,  4 },   /* Aqr */
    { 0x2b2120c0, 0x0ff0,  4 },   /* Aqr */
    { 0x2b210554, 0x1b6c,  3 },   /* Aql */
    { 0x2b2004b0, 0x297c, 16 },   /* Cet */
    { 0x2b0320c0, 0x1f2c,  4 },   /* Aqr */
    { 0x2ae457c6, 0x0d7a, 54 },   /* Mon */
    { 0x2a30fac8, 0x0a8c, 73 },   /* Ser */
    { 0x2a30d41c, 0x10a4, 73 },   /* Ser */
    { 0x2a30ce40, 0x05dc, 47 },   /* Lib */
    { 0x2a3057c6, 0x19e6, 54 },   /* Mon */
    { 0x2a303264, 0x0f3c, 35 },   /* Eri */
    { 0x295e2544, 0x1c5c, 35 },   /* Eri */
    { 0x28aadfd4, 0x1af4, 58 },   /* Oph */
    { 0x28aace40, 0x1194, 47 },   /* Lib */
    { 0x285120c0, 0x2e68,  4 },   /* Aqr */
    { 0x285100a4, 0x08ac, 72 },   /* Sct */
    { 0x2850dfd4, 0x1cd4, 58 },   /* Oph */
    { 0x28505208, 0x1fa4, 54 },   /* Mon */
    { 0x28502544, 0x2238, 35 },   /* Eri */
    { 0x2760972c, 0x0f3c, 26 },   /* Crt */
    { 0x26e94f28, 0x279c, 16 },   /* Cet */
    { 0x2670dfd4, 0x04ec, 71 },   /* Sco */
    { 0x2670c864, 0x1770, 47 },   /* Lib */
    { 0x25f92c00, 0x0780, 11 },   /* Cap */
    { 0x25f91940, 0x0780, 11 },   /* Cap */
    { 0x2580f870, 0x0834, 73 },   /* Ser */
    { 0x2580f168, 0x05dc, 73 },   /* Ser */
    { 0x2508a668, 0x0e10, 28 },   /* Crv */
    { 0x250875a8, 0x2184, 41 },   /* Hya */
    { 0x25086798, 0x0e10, 67 },   /* Pup */
    { 0x25085604, 0x1194, 19 },   /* CMa */
    { 0x25084524, 0x10e0, 46 },   /* Lep */
    { 0x24b8f168, 0x0f3c, 73 },   /* Ser */
    { 0x248d0950, 0x0ff0, 76 },   /* Sgr */
    { 0x236443f8, 0x120c, 46 },   /* Lep */
    { 0x23291940, 0x1a40, 11 },   /* Cap */
    { 0x22b0f780, 0x21c0, 76 },   /* Sgr */
    { 0x22b0e4c0, 0x12c0, 58 },   /* Oph */
    { 0x223875a8, 0x030c, 68 },   /* Pyx */
    { 0x21a2dfd4, 0x0672, 71 },   /* Sco */
    { 0x21487fbc, 0x189c, 41 },   /* Hya */
    { 0x214875a8, 0x0a14, 68 },   /* Pyx */
    { 0x212ae4c0, 0x12c0, 58 },   /* Oph */
    { 0x20d0dc50, 0x0870, 71 },   /* Sco */
    { 0x1fe0b0f4, 0x1770, 41 },   /* Hya */
    { 0x1ef083b8, 0x0564,  1 },   /* Ant */
    { 0x1ef075a8, 0x0e10, 68 },   /* Pyx */
    { 0x1ec21770, 0x1d4c, 36 },   /* For */
    { 0x1eb4891c, 0x48a8, 41 },   /* Hya */
    { 0x1eaadc50, 0x0f3c, 71 },   /* Sco */
    { 0x1e3d4370, 0x2580, 70 },   /* Scl */
    { 0x1e3d2c00, 0x1770, 65 },   /* PsA */
    { 0x1dc483b8, 0x0c6c,  1 },   /* Ant */
    { 0x1d6a4650, 0x0fb4, 22 },   /* Col */
    { 0x1d6a4218, 0x0438,  9 },   /* Cae */
    { 0x1d111df0, 0x0e10, 53 },   /* Mic */
    { 0x1d10f780, 0x2670, 76 },   /* Sgr */
    { 0x1c8483b8, 0x111c,  1 },   /* Ant */
    { 0x1c5cd1c4, 0x0f3c, 49 },   /* Lup */
    { 0x1c5cb0f4, 0x20d0, 14 },   /* Cen */
    { 0x1c20e100, 0x19c8, 71 },   /* Sco */
    { 0x1c204074, 0x05dc,  9 },   /* Cae */
    { 0x1b9483b8, 0x14a0,  1 },   /* Ant */
    { 0x1ab8ac44, 0x2580, 14 },   /* Cen */
    { 0x1ab85c94, 0x1914, 67 },   /* Pup */
    { 0x1ab84650, 0x1644, 22 },   /* Col */
    { 0x19c89ab0, 0x3714, 14 },   /* Cen */
    { 0x19c883b8, 0x16f8,  1 },   /* Ant */
    { 0x19503138, 0x0f3c, 35 },   /* Eri */
    { 0x18f675a8, 0x0e10, 84 },   /* Vel */
    { 0x18d92c00, 0x1c20, 38 },   /* Gru */
    { 0x18d8fac8, 0x12c0, 24 },   /* CrA */
    { 0x18d83c00, 0x0a50,  9 },   /* Cae */
    { 0x17a22a30, 0x11d0, 35 },   /* Eri */
    { 0x178e75a8, 0x2508, 84 },   /* Vel */
    { 0x17714820, 0x2a30, 63 },   /* Phe */
    { 0x17703660, 0x05a0, 40 },   /* Hor */
    { 0x177020d0, 0x1590, 35 },   /* Eri */
    { 0x1680dc50, 0x0a9b, 56 },   /* Nor */
    { 0x1680c738, 0x1518, 49 },   /* Lup */
    { 0x16087080, 0x2a30, 84 },   /* Vel */
    { 0x16085460, 0x1c20, 67 },   /* Pup */
    { 0x160843f8, 0x1068, 64 },   /* Pic */
    { 0x1590300c, 0x0bf4, 40 },   /* Hor */
    { 0x14dd1df0, 0x0e10, 43 },   /* Ind */
    { 0x14dcfd20, 0x20d0, 78 },   /* Tel */
    { 0x14dce6eb, 0x1635,  5 },   /* Ara */
    { 0x14a02a30, 0x11d0, 40 },   /* Hor */
    { 0x14643f48, 0x1518, 64 },   /* Pic */
    { 0x13b0d7a0, 0x0f4b, 56 },   /* Nor */
    { 0x139c19c8, 0x1068, 35 },   /* Eri */
    { 0x1338396c, 0x05dc, 32 },   /* Dor */
    { 0x13382580, 0x13ec, 40 },   /* Hor */
    { 0x12c11df0, 0x1770, 43 },   /* Ind */
    { 0x12665460, 0x1e78, 12 },   /* Car */
    { 0x124835e8, 0x0960, 32 },   /* Dor */
    { 0x124821fc, 0x13ec, 40 },   /* Hor */
    { 0x120c1644, 0x0bb8, 35 },   /* Eri */
    { 0x11943f48, 0x1770, 64 },   /* Pic */
    { 0x115856b8, 0x201c, 12 },   /* Car */
    { 0x11443138, 0x0708, 69 },   /* Ret */
    { 0x111c12c0, 0x0f3c, 35 },   /* Eri */
    { 0x10e0d3a4, 0x1347, 56 },   /* Nor */
    { 0x10e03840, 0x0e10, 32 },   /* Dor */
    { 0x10e01e78, 0x12c0, 40 },   /* Hor */
    { 0x10a456b8, 0x2580, 12 },   /* Car */
    { 0x1068cc60, 0x0b40, 18 },   /* Cir */
    { 0x1068b478, 0x17e8, 14 },   /* Cen */
    { 0x1068a668, 0x0e10, 27 },   /* Cru */
    { 0x10684650, 0x1518, 64 },   /* Pic */
    { 0x0fb45b68, 0x42cc, 12 },   /* Car */
    { 0x0fb43138, 0x0bb8, 69 },   /* Ret */
    { 0x0f793560, 0x12c0, 81 },   /* Tuc */
    { 0x0f78f618, 0x27d8, 60 },   /* Pav */
    { 0x0f3c3cf0, 0x1068, 32 },   /* Dor */
    { 0x0f3c2d00, 0x0ff0, 69 },   /* Ret */
    { 0x0f004d58, 0x12c0, 64 },   /* Pic */
    { 0x0ec53560, 0x2ee0, 81 },   /* Tuc */
    { 0x0ec412c0, 0x0bb8, 42 },   /* Hyi */
    { 0x0e882d00, 0x1374, 69 },   /* Ret */
    { 0x0e10f618, 0x35e8, 60 },   /* Pav */
    { 0x0e10d7a0, 0x0f4b, 79 },   /* TrA */
    { 0x0d98d548, 0x13ec, 79 },   /* TrA */
    { 0x0d984074, 0x13ec, 32 },   /* Dor */
    { 0x0c62d1c4, 0x19c8, 79 },   /* TrA */
    { 0x0c30bdd8, 0x13ec, 18 },   /* Cir */
    { 0x0c309e34, 0x1fa4, 55 },   /* Mus */
    { 0x0c305c94, 0x2274, 86 },   /* Vol */
    { 0x0c304074, 0x1c20, 32 },   /* Dor */
    { 0x0bb8d1c4, 0x1af4, 79 },   /* TrA */
    { 0x0bb89e34, 0x21fc, 55 },   /* Mus */
    { 0x0a8d2c00, 0x1c20, 43 },   /* Ind */
    { 0x0a8cef10, 0x0e10,  2 },   /* Aps */
    { 0x0a8ccf6c, 0x1fa4, 79 },   /* TrA */
    { 0x0a8c12c0, 0x2db4, 42 },   /* Hyi */
    { 0x0960c030, 0x3cf0,  2 },   /* Aps */
    { 0x09604074, 0x1c20, 52 },   /* Men */
    { 0x0708fd20, 0x5460, 57 },   /* Oct */
    { 0x07086bd0, 0x5460, 17 },   /* Cha */
    { 0x07083138, 0x3a98, 52 },   /* Men */
    { 0x07080000, 0x0a8c, 42 },   /* Hyi */
    { 0x06900000, 0x3138, 42 },   /* Hyi */
    { 0x0384f744, 0x8b74, 57 },   /* Oct */
    { 0x03846bd0, 0x8b74, 57 },   /* Oct */
    { 0x0258a8c0, 0xa8c0, 57 },   /* Oct */
    { 0x02580000, 0xa8c0, 57 } }; /* Oct */

/* Each of the above represents a horizontal constellation boundary,
at south polar distance (in arcminutes) SPD(index),  minimum RA of
RA(index) (in seconds of RA),  extending by d_ra seconds,  with the
indicated constellation to the south.  As the macros indicate,  RA
is stored in the lower 17 bits of 'ra_spd',  and SPD in the next 14
bits (the uppermost bit is unused).  */

#define SPD(i)   (bounds[i].ra_spd >> 17)
#define RA(i)  (bounds[i].ra_spd & 0x1ffff)

#ifdef TEST_MAIN
   #include <stdio.h>
   #include <stdlib.h>
#endif

int DLL_FUNC constell_from_ra_dec( const double ra_degrees_1875,
                                   const double dec_degrees_1875,
                                   char DLLPTR *constell_name)
{
   const int32_t ra = ((int32_t)(ra_degrees_1875 * 240.) + 864000) % 86400;
   const int16_t spd = (int16_t)( (dec_degrees_1875 + 90.) * 60.);
   int idx = -1, step = 512;
   int rval = -1;
   const int n_bounds = (int)( sizeof( bounds) / sizeof( bounds[0]));      /* = 324 */

   while( step >>= 1)
      if( idx + step < n_bounds && spd < SPD( idx + step))
         idx += step;
                  /* 'idx' = index of first entry north of dec.  May be -1 */
                  /* if dec > 88 and there are no entries north of us.  */
   while( idx >= 0 && rval == -1)
      {
      const int32_t min_ra = RA( idx);
      const int32_t max_ra = min_ra + (int32_t)bounds[idx].d_ra;

#ifdef TEST_MAIN
      printf( "%4d %11.8f %11.8f %+011.7f %.3s\n", idx, (double)min_ra / 3600.,
                           (double)max_ra / 3600., (double)SPD( idx) / 60. - 90.,
                           constell_names + 3 * bounds[idx].constell_idx);
#endif
      if( (ra >= min_ra && ra < max_ra) || (ra + 86400 >= min_ra && ra + 86400 < max_ra))
         rval = bounds[idx].constell_idx;
      idx--;
      }
   if( rval == -1)    /* Didn't hit anything;  must be in UMi */
      rval = 83;
   if( constell_name)
      {
      memcpy( constell_name, constell_names + 3 * rval, 3);
      constell_name[3] = '\0';
      }
   return( rval);
}

#ifdef TEST_MAIN

/* Run this with an RA/dec in decimal hours and degrees.  The corresponding
constellation and some of the logic required to determine it will be shown. */

int main( const int argc, const char **argv)
{
   if( argc == 3)
      {
      char constell[4];
      int rval = constell_from_ra_dec( atof( argv[1]) * 15, atof( argv[2]), constell);

      printf( "Memory used: %d (%d bytes/entry)\n", (int)sizeof( bounds),
                              (int)sizeof( bounds[0]));
      printf( "rval %d = '%s'\n", rval, constell);
      }
   else
      printf( "Test/debugging code to determine constellation for a given point\n"
              "Usage : conbound (RA in decimal hours) (dec in decimal degrees)\n"
              "Coordinates must be in B1875 epoch\n");
   return( 0);
}
#endif
