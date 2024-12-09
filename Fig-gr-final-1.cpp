//"From experimental clues to theoretical modeling: Evolution associated with the membrane-takeover at an early stage of life"
//by Wentao Ma*, Chunwu Yu
//C source codes for the simulation program  --- The case is correponding to Fig.2a in the article

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <conio.h>

/******* RANDOM NUMBER GENERATOR ********/
#define A1 471
#define B1 1586
#define C1 6988
#define D1 9689
#define M 16383
#define RIMAX 2147483648.0        // = 2^31 
#define RandomInteger (++nd, ra[nd & M] = ra[(nd-A1) & M] ^ ra[(nd-B1) & M] ^ ra[(nd-C1) & M] ^ ra[(nd-D1) & M])
static long ra[M + 1], nd;
void seed(long seed);  // Random number initialization 
long randl(long);      // Random number between 0 and parameter 
double randd(void);    // Random double between 0 and 1         
/****************************************/

#define LEN sizeof(struct rna)
#define C 2
#define G 3
#define A 1
#define U 4
#define STEPNUM 2000000    // Total time steps of Monte Carlo simulation
#define STAREC 0           // The step to start record
#define RECINT 10000       // The interval steps of recording
#define MAX_RNA_LENGTH 100  // Maximum RNA length allowed in the simulation
#define MAXRNASINGRID 1000  // Maximum number of RNA allowed in a single grid room

#define NRSEQ G,C,A,C,G,U,A   // The presumed characteristic sequence of a nucleotide-synthetase ribozyme (NR)
#define GRSEQ U,U,G,A,G,C,G   // The presumed characteristic sequence of a glycerophosphate-synthetase ribozyme (GR)
#define NPRSEQ U,C,A,C,G,A,G  // The presumed characteristic sequence of a nucleotide-precursor-synthetase ribozyme (NPR)
#define CONTRSEQ1 G,G,C,U,A,C,U  // The presumed characteristic sequence of control-1
#define CONTRSEQ2 U,G,A,C,C,A,U  // The presumed characteristic sequence of control-2

#define INOCUSTEP0 1000       // The step of the first inoculation
#define INOCUSTEP1 10000      // The step of the second inoculation

#define MECH_CHANGE_STEP 0  // The step of changing relevant mechanism

#define INOCU_E_CELL_NUM 1     // The number of inoculated empty protocells 
#define INOCU_CELL_NUM 10      // The number of inoculated RNA-based protocells
#define INOCU_SEQ_NUM 1        // The number of inoculated functional RNA molecules in a protocell
#define INOCUSEQ1 GRSEQ        
#define INOCUSEQ2 NRSEQ
#define INOCUSEQ3 CONTRSEQ1
#define INOCUSEQ4 CONTRSEQ2
#define INOCUSEQ5 NPRSEQ

#define SD 7           // The random seed

#define N 30          // The system surface is defined as an N×N grid
#define ROOMNUM (N*N) // The room number in the grid
#define TNPPB 80000   // Total nucleotide precursors' precursor (quotients in measurement of nucleotides) introduced in the beginning
#define TFB   50000   // Total fatty acids introduced in the beginning
#define TGPB  50000   //Total glycerophosphate precursors (quotients in measurement of glycerophosphate) introduced in the beginning

#define PNF   0.005   // Probability of a nucleotide forming from its precursor (not catalyzed)
#define PNFR  0.2     // Probability of a nucleotide forming from its precursor catalyzed by the nucleotide synthetase ribozyme
#define NR_TURN 1     // The catalytic turns of an NR in a time step
#define PNPF  0.002   // Probability of a nucleotide-precursor forming from its precursor (not catalyzed)
#define PNPFR 0.0     // Probability of a nucleotide-precursor forming from its precursor catalyzed by NPR
#define NPR_TURN 1    // The catalytic turns of an NPR in a time step
#define PGF   0.002   // Probability of a glycerophosphate forming from its precursor (not catalyzed)
#define PGFR  0.9     // Probability of a glycerophosphate forming from its precursor catalyzed by GR
#define GR_TURN 1     // The catalytic turns of a GR in a time step
#define PPF   0.02    // Probability of a phospholipid forming from a glycerophosphate and two fatty acids in situ on the membrane

#define PND   0.02   // Probability of a nucleotide decaying into its precursor
#define PNDE  0.001  // Probability of a nucleotide residue decaying at RNA's chain end
#define PNPD  0.005  // Probability of a nucleotide precursor decaying into its precursor
#define PGD   0.1    // Probability of a glycerophosphate decaying into its precursor
#define PPD   0.1    // Probability of a phospholipid decay (into fatty acid and glycerophosphate) out of the membrane
#define PPDM  0.01   // Probability of a phospholipid decay within the membrane

#define PRL 0.000001 // Probability of the random ligation of nucleotides and RNA
#define PBB 0.00001  // Probability of a phosphodiester bond breaking within an RNA chain

#define PAT  0.9     // Probability of an RNA template attracting a substrate (by base-pairing)
#define PFP  0.001   // Probability of the false base-pairing (related to RNA's replicating fidelity)
#define PTL  0.02    // Probability of the template-directed ligation of RNA
#define PSP  0.5    // Probability of the separation of a base pair

#define PFLM 0.002  // Probability of a fatty acid leaving the membrane
#define PPLM 0.0001 // Probability of a phospholipid leaving the membrane
#define PFJM 0.9    // Probability of a fatty acid joining the membrane
#define PPJM 0.9    // Probability of a phospholipid joining the membrane

#define PNP 0.00005 // Probability of a nucleotide permeating through the membrane
#define PNPP  0.05  // Probability of a nucleotide precursor permeating through the membrane
#define PNPPP 0.5   // Probability of a nucleotide-precursor's precursor permeating through the membrane
#define PGPP  0.9   // Probability of a glycerophosphate precursor permeating through the membrane

#define PMF 0.1     // Probability of a membrane forming  
#define PCB 0.0001  // Probability of a protocell breaking   
#define PCD 0.1     // Probability of a protocell dividing
#define PCF 0.001   // Probability of two adjacent protocells fusing with each other

#define PMC 0.1     // Probability of the movement of protocells
#define PMV 0.9     // Probability of the movement of a nucleotide/fatty acid/phospholipid (or relevant precursors)

#define FDE  1.0    // Factor of the Donnan's equilibrium effect
#define FDO  10     // Factor of molecular degradation outside protocells
#define FPL  5      // Factor of phospholipids' influence on amphiphiles leaving the membrane
#define FPP  20     // Factor of phospholipids on permeability (for nucleotides or their precursors)
#define FPPW 3      // The Weak version of FPP (for nucleotide-precursor's precursors or glycerophosphate precursors)
#define LAM  200.0  // The lower limit number of amphiphiles to form a protocell membrane

#define RMRW  (pow(p->length1+p->length2,1/2.0))  // The relationship between the movement of an RNA molecule and its molecular weight

void avail_xy_init(void); // Initialization for xy_choose
void xy_choose(void);     // Picks a room at random 
void fresh_rna(int h);    // Updating a unit for the next time step
int findseq(char subseq[], int subseqlength, struct rna* p); //find a specific subsequence in an RNA 
void cell_fusing(int ay, int ax); // Fusion between two adjacent protocells
void oneway_outer_moving(int ay, int by, int ax, int bx); // The oneway movement of outer molecules when a protocell moves or divides
void twoway_outer_moving(int ay, int by, int cy, int ax, int bx, int cx); // The twoway movement of outer molecules when a protocell moves or divides
void threeway_outer_moving(int ay, int by, int cy, int dy, int ax, int bx, int cx, int dx); // The threeway movement of outer molecules when a protocell moves or divides
void cell_moving(int ay, int ax);    // Movement of protocells
void cell_dividing(int ay, int ax);  // Division of protocells
void inits(void);         // Initialization of the system
void inocu_cell(void);    // Inoculation of protocells
void inoculate1(void);    // Inoculation of RNA
void unit_action(void);   // Action of units (molecules and protocells) in the system
void record(void);        // Data recording at every interval step (RECINT) 
void freepool(void);      // Memory releasing

struct rna                // Nucleotide or RNA
{
	char information[2][MAX_RNA_LENGTH];
	int length1;
	int length2;
	int nick[MAX_RNA_LENGTH];
	struct rna* next;
	struct rna* prior;
};
struct rna* room_head[2][N][N]; // Nucleotides or RNA in grid rooms
struct rna* p, * p1, * p2, * p3, * p4, * ps, * ps1, * ps2;
struct rna* rna_arr[MAXRNASINGRID];

static char nrseq[50] = { NRSEQ };
static char grseq[50] = { GRSEQ };
static char nprseq[50] = { NPRSEQ };
static char contrseq1[50] = { CONTRSEQ1 };
static char contrseq2[50] = { CONTRSEQ2 };
static char inocuseq1[50] = { INOCUSEQ1 };
static char inocuseq2[50] = { INOCUSEQ2 };
static char inocuseq3[50] = { INOCUSEQ3 };
static char inocuseq4[50] = { INOCUSEQ4 };
static char inocuseq5[50] = { INOCUSEQ5 };
static int np_arr[2][N][N];  // Nucleotide precursors in grid rooms
static int npp_arr[2][N][N]; // Nucleotide-precursor's precursors in grid rooms
static int fa_arr[2][N][N];  // Fatty acids in grid rooms

static int gp_arr[2][N][N];  // Glycerophosphate precursors in grid rooms
static int g_arr[2][N][N];   // Glycerophosphates in grid rooms
static int phl_arr[2][N][N]; // Phospholipids in grid rooms

static int m_arr[N][N];      // Fatty acids and their residues (in phospholipids) on the membrane of protocells in grid rooms
static int phl_on_m_arr[N][N];   // Phospholipids on the membrane of protocells in grid rooms

int flagcmov[N][N];    // Flags labeling the movement of protocells in grid rooms
int flagcdiv[N][N];    // Flags labeling the devision of protocells in grid rooms

int x, y;                // The coordinate of rooms in the grid 
long i;                  // Cycle variable for Monte Carlo steps
int g = 0;				 // Recording times
int over_max_len = 0;    // Record the times that RNA may exceed the  maximum length setting 
long cell_div_times = 0;
long available;
long availabl[ROOMNUM];
int k, np_bef, npp_bef, fa_bef, gp_bef, g_bef, phl_bef, m_bef, m_fcd_bef, m_phl_bef, randcase;
int nrlength, nprlength, contrlength1, contrlength2, grlength;
int inoculength1, inoculength2, inoculength3, inoculength4, inoculength5;
int flag1, flag2, flag3, flagnr, flagnpr, flaggr, flagctrl, flagntsyn, flagnpsyn, flaggsyn;

float gr_num[(STEPNUM - STAREC) / RECINT + 1];      // For recording the number of GR
float nr_num[(STEPNUM - STAREC) / RECINT + 1];      // For recording the number of NR
float npr_num[(STEPNUM - STAREC) / RECINT + 1];     // For recording the number of NPR
float contr1_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of control-1
float contr2_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of control-2
float gr_incell_num[(STEPNUM - STAREC) / RECINT + 1];   // For recording the number of GR inside protocells 
float nr_incell_num[(STEPNUM - STAREC) / RECINT + 1];   // For recording the number of NR inside protocells 
float npr_incell_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of NPR inside protocells 
float contr1_incell_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of control-1 inside protocells 
float contr2_incell_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of control-2 inside protocells

float total_nt_mat[(STEPNUM - STAREC) / RECINT + 1];  // For recording the total nucleotide materials 
float total_fa_mat[(STEPNUM - STAREC) / RECINT + 1];  // For recording the total fatty acid materials 
float total_g_mat[(STEPNUM - STAREC) / RECINT + 1];   // For recording the total glycerophosphate materials 

float rna_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of nucleotides and RNA
float np_num[(STEPNUM - STAREC) / RECINT + 1];   // For recording the number of nucleotide precursors
float npp_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of nucleotide-precursor's precursors 
float fa_num[(STEPNUM - STAREC) / RECINT + 1];   // For recording the number of fatty acids

float gp_num[(STEPNUM - STAREC) / RECINT + 1];   // For recording the number of glycerophosphate precursors
float g_num[(STEPNUM - STAREC) / RECINT + 1];    // For recording the number of glycerol-phosphates
float phl_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of phospholipids
float phl_on_m_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of phospholipids on the membrane

float cell_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of protocells
float cell_0_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of protocells without functional RNA species
float cell_gr_num[(STEPNUM - STAREC) / RECINT + 1]; // For recording the number of protocells containing GR
float cell_nr_num[(STEPNUM - STAREC) / RECINT + 1]; // For recording the number of protocells containing NR
float cell_npr_num[(STEPNUM - STAREC) / RECINT + 1]; // For recording the number of protocells containing NPR
float cell_nrgr_num[(STEPNUM - STAREC) / RECINT + 1]; // For recording the number of protocells containing both NR and GR 
float cell_nrnpr_num[(STEPNUM - STAREC) / RECINT + 1]; // For recording the number of protocells containing both NR and NPR 
float cell_grnpr_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of protocells containing both GR and NPR 
float cell_nrgrnpr_num[(STEPNUM - STAREC) / RECINT + 1]; // For recording the number of protocells containing all NR, GR and NPR 

float cell_ctrl_num[(STEPNUM - STAREC) / RECINT + 1];  // For recording the number of protocells containing control-1

/////////////////////////////////////////////////////////////////////////
void seed(long seed)   // Random number initialization    
{
	int a;
	if (seed < 0) { puts("SEED error."); exit(1); }
	ra[0] = (long)fmod(16807.0 * (double)seed, 2147483647.0);
	for (a = 1; a <= M; a++)
	{
		ra[a] = (long)fmod(16807.0 * (double)ra[a - 1], 2147483647.0);
	}
}

/////////////////////////////////////////////////////////////////////////
long randl(long num)      // Random integer number between 0 and num-1 
{
	return(RandomInteger % num);
}

/////////////////////////////////////////////////////////////////////////
double randd(void)        // Random real number between 0 and 1 
{
	return((double)RandomInteger / RIMAX);
}

/////////////////////////////////////////////////////////////////////////
void avail_xy_init(void)   // Initialization for xy_choose
{
	int j;
	for (j = 0; j < ROOMNUM; j++)
	{
		availabl[j] = j + 1;
	}
	available = ROOMNUM;
}

/////////////////////////////////////////////////////////////////////////
void xy_choose(void)       // Picks a room at random
{
	long rl, s;
	rl = randl(available);
	s = availabl[rl];
	x = (s - 1) % N;
	y = (s - 1) / N;
	availabl[rl] = availabl[available - 1];
	available--;
}

/////////////////////////////////////////////////////////////////////////
void fresh_rna(int h)     // Updating a nucleotide or RNA for the next time step
{
	p1 = p->prior;
	p2 = p->next;
	p3 = room_head[!h][y][x]->next;
	room_head[!h][y][x]->next = p;
	p->next = p3;
	p->prior = room_head[!h][y][x];
	if (p3 != room_head[!h][y][x])p3->prior = p;
	p1->next = p2;
	if (p2 != room_head[h][y][x])p2->prior = p1;
	p = p1;
}

/////////////////////////////////////////////////////////////////////////
void rna_shuffle(int h) //Random order-arrangement of nodes in the linked list for nucleotides and RNA in a grid room 
{
	int a, b, len;
	for (a = 0, p = room_head[h][y][x]->next; p != room_head[h][y][x]; p = p->next)
	{
		if (a == MAXRNASINGRID) { printf("Too many RNA molecules in a single grid room"); exit(0); }
		rna_arr[a] = p;
		a++;
	}
	len = a;
	rna_arr[a] = NULL;
	for (a = len - 1; a >= 1; a--)
	{
		p = rna_arr[a];
		b = RandomInteger % (a);
		rna_arr[a] = rna_arr[b];
		rna_arr[b] = p;
	}
	p = room_head[h][y][x];
	for (a = 0; a < len; a++)
	{
		p->next = rna_arr[a];
		rna_arr[a]->prior = p;
		p = p->next;
	}
	p->next = room_head[h][y][x];

}
/////////////////////////////////////////////////////////////////////////
int findseq(char subseq[], int subseqlength, struct rna* p)  // Find a specific subsequence in an RNA
{
	int a, b, flag1, flag2;

	flag1 = 0;
	if (p->length1 >= subseqlength)
	{
		for (b = 0; p->length1 - subseqlength - b >= 0; b++)
		{
			flag2 = 0;
			for (a = 0; a < subseqlength; a++)
			{
				if (p->information[0][b + a] == subseq[a])continue;
				else { flag2 = 1; break; }
			}
			if (flag2 == 0)break;
		}
		if (flag2 == 1)flag1 = 1;
	}
	else flag1 = 1;

	if (flag1 == 0)return(0);   // Yes, the RNA contains the subsequence
	else return(1);   // No, the RNA does not contain the subsequence
}

/////////////////////////////////////////////////////////////////////////
void cell_fusing(int ay, int ax)
{
	fa_arr[0][(N + y + ay) % N][(N + x + ax) % N] += fa_arr[0][y][x]; // Inner fatty acids moving 
	fa_arr[0][y][x] = 0;

	g_arr[0][(N + y + ay) % N][(N + x + ax) % N] += g_arr[0][y][x]; // Inner glycerophosphates moving 
	g_arr[0][y][x] = 0;

	gp_arr[0][(N + y + ay) % N][(N + x + ax) % N] += gp_arr[0][y][x]; // Inner glycerophosphate precursors moving 
	gp_arr[0][y][x] = 0;

	phl_arr[0][(N + y + ay) % N][(N + x + ax) % N] += phl_arr[0][y][x]; // Inner phospholipids moving 
	phl_arr[0][y][x] = 0;

	np_arr[0][(N + y + ay) % N][(N + x + ax) % N] += np_arr[0][y][x];  // Inner nucleotide precursors moving 
	np_arr[0][y][x] = 0;

	npp_arr[0][(N + y + ay) % N][(N + x + ax) % N] += npp_arr[0][y][x]; // Inner nucleotide-precursor's precursors moving 
	npp_arr[0][y][x] = 0;

	if (room_head[0][y][x]->next != room_head[0][y][x])  // Inner nucleotides and RNA moving
	{
		p = room_head[0][y][x];
		do {
			p = p->next;
		} while (p->next != room_head[0][y][x]); // Find the end of the linked list
		p->next = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next;
		if (p->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]) p->next->prior = p;
		room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = room_head[0][y][x]->next;
		(room_head[0][y][x]->next)->prior = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
		room_head[0][y][x]->next = room_head[0][y][x];
	}

	m_arr[(N + y + ay) % N][(N + x + ax) % N] += m_arr[y][x]; // Membrane fusing 
	phl_on_m_arr[(N + y + ay) % N][(N + x + ax) % N] += phl_on_m_arr[y][x];
	m_arr[y][x] = 0;
	phl_on_m_arr[y][x] = 0;
}

/////////////////////////////////////////////////////////////////////////
void oneway_outer_moving(int ay, int by, int ax, int bx)
{
	fa_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N] += fa_arr[0][(N + y + ay) % N][(N + x + ax) % N]; // Outer fatty acids moving
	fa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;

	g_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N] += g_arr[0][(N + y + ay) % N][(N + x + ax) % N]; // Outer glycerophosphates moving
	g_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;

	gp_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N] += gp_arr[0][(N + y + ay) % N][(N + x + ax) % N]; // Outer glycerol-phosphate precusors moving
	gp_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;

	phl_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N] += phl_arr[0][(N + y + ay) % N][(N + x + ax) % N]; // Outer phospholipids moving
	phl_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;

	np_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N] += np_arr[0][(N + y + ay) % N][(N + x + ax) % N];  // Outer nucleotide precursors moving
	np_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;

	npp_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N] += npp_arr[0][(N + y + ay) % N][(N + x + ax) % N];  // Outer nucleotide-precursor's precursors moving
	npp_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;

	if (room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]) // Outer nucleotides and RNA moving
	{
		p = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
		do {
			p = p->next;
		} while (p->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]);
		p->next = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next;
		if (room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next != room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N])(room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next)->prior = p;
		room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next;
		(room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)->prior = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N];
		room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
	}
}

/////////////////////////////////////////////////////////////////////////
void twoway_outer_moving(int ay, int by, int cy, int ax, int bx, int cx)
{
	fa_bef = fa_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer fatty acids moving
	fa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < fa_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  fa_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  fa_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		default: printf("two way outer fa moving error");
		}
	}

	g_bef = g_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer glycerophosphates moving
	g_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < g_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  g_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  g_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		default: printf("two way outer g moving error");
		}
	}

	gp_bef = gp_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer glycerophosphate precusors moving
	gp_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < gp_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  gp_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  gp_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		default: printf("two way outer gp moving error");
		}
	}

	phl_bef = phl_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer phospholipids moving
	phl_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < phl_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  phl_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  phl_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		default: printf("two way outer pl moving error");
		}
	}

	np_bef = np_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer nucleotide precursors moving
	np_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < np_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  np_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  np_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		default: printf("two way outer np moving error");
		}
	}

	npp_bef = npp_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer nucleotide-precursor's precursors moving
	npp_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < npp_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  npp_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  npp_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		default: printf("two way outer npp moving error");
		}
	}

	for (p = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next; p != room_head[0][(N + y + ay) % N][(N + x + ax) % N]; p = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)
	{	                                                     // Outer nucleotides and RNA moving
		room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = p->next;
		if (p->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]) (p->next)->prior = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
		randcase = randl(2);
		switch (randcase)
		{
		case 0:
			p->next = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next;
			if (room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next != room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]) (room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next)->prior = p;
			room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next = p;
			p->prior = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N];
			break;
		case 1:
			p->next = room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next;
			if (room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next != room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]) (room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next)->prior = p;
			room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next = p;
			p->prior = room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N];
			break;
		default: printf("twoway outer RNA moving error");
		}
	}
}

/////////////////////////////////////////////////////////////////////////
void threeway_outer_moving(int ay, int by, int cy, int dy, int ax, int bx, int cx, int dx)
{
	fa_bef = fa_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer fatty acids moving
	fa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < fa_bef; k++)
	{
		randcase = randl(3);
		switch (randcase)
		{
		case 0:  fa_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  fa_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		case 2:  fa_arr[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]++; break;
		default: printf("three way outer fa moving error");
		}
	}

	g_bef = g_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer glycerophosphates moving
	g_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < g_bef; k++)
	{
		randcase = randl(3);
		switch (randcase)
		{
		case 0:  g_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  g_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		case 2:  g_arr[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]++; break;
		default: printf("three way outer g moving error");
		}
	}

	gp_bef = gp_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer glycerophosphate precursors moving
	gp_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < gp_bef; k++)
	{
		randcase = randl(3);
		switch (randcase)
		{
		case 0:  gp_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  gp_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		case 2:  gp_arr[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]++; break;
		default: printf("three way outer gp moving error");
		}
	}

	phl_bef = phl_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer phospholipids moving
	phl_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < phl_bef; k++)
	{
		randcase = randl(3);
		switch (randcase)
		{
		case 0:  phl_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  phl_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		case 2:  phl_arr[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]++; break;
		default: printf("three way outer pl moving error");
		}
	}

	np_bef = np_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer nucleotide precursors moving
	np_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < np_bef; k++)
	{
		randcase = randl(3);
		switch (randcase)
		{
		case 0:  np_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  np_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		case 2:  np_arr[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]++; break;
		default: printf("three way outer np moving error");
		}
	}

	npp_bef = npp_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer nucleotide-precursor's precursors moving
	npp_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
	for (k = 0; k < npp_bef; k++)
	{
		randcase = randl(3);
		switch (randcase)
		{
		case 0:  npp_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
		case 1:  npp_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
		case 2:  npp_arr[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]++; break;
		default: printf("three way outer npp moving error");
		}
	}

	for (p = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next; p != room_head[0][(N + y + ay) % N][(N + x + ax) % N]; p = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)
	{	                                                    // Outer nucleotides and RNA moving
		room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = p->next;
		if (p->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]) (p->next)->prior = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
		randcase = randl(3);
		switch (randcase)
		{
		case 0:
			p->next = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next;
			if (room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next != room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]) (room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next)->prior = p;
			room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next = p;
			p->prior = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N];
			break;
		case 1:
			p->next = room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next;
			if (room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next != room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]) (room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next)->prior = p;
			room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next = p;
			p->prior = room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N];
			break;
		case 2:
			p->next = room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next;
			if (room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next != room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]) (room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next)->prior = p;
			room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next = p;
			p->prior = room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N];
			break;
		default: printf("threeway outer RNA moving error");
		}
	}
}

/////////////////////////////////////////////////////////////////////////
void cell_moving(int ay, int ax)
{
	fa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = fa_arr[0][y][x];  // Inner fatty acids moving 
	fa_arr[0][y][x] = 0;

	g_arr[0][(N + y + ay) % N][(N + x + ax) % N] = g_arr[0][y][x];   // Inner glycerophosphates moving 
	g_arr[0][y][x] = 0;

	gp_arr[0][(N + y + ay) % N][(N + x + ax) % N] = gp_arr[0][y][x]; // Inner glycerophosphate precursors moving 
	gp_arr[0][y][x] = 0;

	phl_arr[0][(N + y + ay) % N][(N + x + ax) % N] = phl_arr[0][y][x]; // Inner phospholipids moving 
	phl_arr[0][y][x] = 0;

	np_arr[0][(N + y + ay) % N][(N + x + ax) % N] = np_arr[0][y][x];  // Inner nucleotide precursors moving 
	np_arr[0][y][x] = 0;

	npp_arr[0][(N + y + ay) % N][(N + x + ax) % N] = npp_arr[0][y][x]; // Inner nucleotide-precursor's precursors moving 
	npp_arr[0][y][x] = 0;

	if (room_head[0][y][x]->next != room_head[0][y][x])  // Inner nucleotides and RNA moving
	{
		p = room_head[0][y][x];
		do {
			p = p->next;
		} while (p->next != room_head[0][y][x]);
		p->next = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
		room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = room_head[0][y][x]->next;
		(room_head[0][y][x]->next)->prior = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
		room_head[0][y][x]->next = room_head[0][y][x];
	}

	m_arr[(N + y + ay) % N][(N + x + ax) % N] = m_arr[y][x]; // Membrane moving 
	phl_on_m_arr[(N + y + ay) % N][(N + x + ax) % N] = phl_on_m_arr[y][x];
	m_arr[y][x] = 0;
	phl_on_m_arr[y][x] = 0;
}

/////////////////////////////////////////////////////////////////////////
void cell_dividing(int ay, int ax)
{
	fa_bef = fa_arr[0][y][x];  // Inner fatty acids distributing
	for (k = 0; k < fa_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  // Staying in the old cell
			break;
		case 1:  // Moving to the new cell 
			fa_arr[0][(N + y + ay) % N][(N + x + ax) % N]++; fa_arr[0][y][x]--; break;
		default:
			printf("inner fa moving error in cell division");
		}
	}

	g_bef = g_arr[0][y][x];  // Inner glycerophosphates distributing
	for (k = 0; k < g_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  // Staying in the old cell
			break;
		case 1:  // Moving to the new cell 
			g_arr[0][(N + y + ay) % N][(N + x + ax) % N]++; g_arr[0][y][x]--; break;
		default:
			printf("inner g moving error in cell division");
		}
	}

	gp_bef = gp_arr[0][y][x];  // Inner glycerophosphate precursors distributing
	for (k = 0; k < gp_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  // Staying in the old cell
			break;
		case 1:  // Moving to the new cell 
			gp_arr[0][(N + y + ay) % N][(N + x + ax) % N]++; gp_arr[0][y][x]--; break;
		default:
			printf("inner gp moving error in cell division");
		}
	}

	phl_bef = phl_arr[0][y][x];  // Inner phospholipids distributing
	for (k = 0; k < phl_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  // Staying in the old cell
			break;
		case 1:  // Moving to the new cell 
			phl_arr[0][(N + y + ay) % N][(N + x + ax) % N]++; phl_arr[0][y][x]--; break;
		default:
			printf("inner pl moving error in cell division");
		}
	}

	np_bef = np_arr[0][y][x];  // Inner nucleotide precursors distributing
	for (k = 0; k < np_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  // Staying in the old cell
			break;
		case 1:  // Moving to the new cell 
			np_arr[0][(N + y + ay) % N][(N + x + ax) % N]++; np_arr[0][y][x]--; break;
		default:
			printf("inner np moving error in cell division");
		}
	}

	npp_bef = npp_arr[0][y][x];  // Inner nucleotide-precursor's precursors distributing
	for (k = 0; k < npp_bef; k++)
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:  // Staying in the old cell
			break;
		case 1:  // Moving to the new cell 
			npp_arr[0][(N + y + ay) % N][(N + x + ax) % N]++; npp_arr[0][y][x]--; break;
		default:
			printf("inner npp moving error in cell division");
		}
	}

	for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next) // Inner nucleotides and RNA distributing
	{
		randcase = randl(2);
		switch (randcase)
		{
		case 0:   // Staying in the old cell
			break;
		case 1:   // Moving to the new cell
			(p->prior)->next = p->next;
			if (p->next != room_head[0][y][x]) (p->next)->prior = p->prior;
			p1 = p;
			p = p->prior;
			p1->next = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next;
			if (room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]) (room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)->prior = p1;
			room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = p1;
			p1->prior = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
			break;
		default: printf("inner RNA moving error");
		}
	}

	// Membrane dividing
	m_fcd_bef = m_arr[y][x] - 2 * phl_on_m_arr[y][x];
	m_phl_bef = phl_on_m_arr[y][x];

	m_arr[y][x] = 0;
	phl_on_m_arr[y][x] = 0;
	m_arr[(N + y + ay) % N][(N + x + ax) % N] = 0;
	phl_on_m_arr[(N + y + ay) % N][(N + x + ax) % N] = 0;

	flag1 = 0; flag2 = 0;  // Labels for an offspring protocell's membrane reaching LAM
	while (m_fcd_bef + m_phl_bef)
	{
		if (m_fcd_bef > 0)
		{
			m_fcd_bef--;
			if (flag1 == 1 && flag2 == 0)
			{
				m_arr[(N + y + ay) % N][(N + x + ax) % N]++;
				if (m_arr[(N + y + ay) % N][(N + x + ax) % N] >= LAM) flag2 = 1;
			}
			else if (flag1 == 0 && flag2 == 1)
			{
				m_arr[y][x]++;
				if (m_arr[y][x] >= LAM) flag1 = 1;
			}
			else
			{
				if (randd() < 0.5)
				{
					m_arr[y][x]++;
					if (m_arr[y][x] >= LAM) flag1 = 1;
				}
				else
				{
					m_arr[(N + y + ay) % N][(N + x + ax) % N]++;
					if (m_arr[(N + y + ay) % N][(N + x + ax) % N] >= LAM) flag2 = 1;
				}
			}
		}
		if (m_phl_bef > 0)
		{
			m_phl_bef--;
			if (flag1 == 1 && flag2 == 0)
			{
				m_arr[(N + y + ay) % N][(N + x + ax) % N] += 2;
				phl_on_m_arr[(N + y + ay) % N][(N + x + ax) % N]++;
				if (m_arr[(N + y + ay) % N][(N + x + ax) % N] >= LAM) flag2 = 1;
			}
			else if (flag1 == 0 && flag2 == 1)
			{
				m_arr[y][x] += 2;
				phl_on_m_arr[y][x]++;
				if (m_arr[y][x] >= LAM) flag1 = 1;
			}
			else
			{
				if (randd() < 0.5)
				{
					m_arr[y][x] += 2; phl_on_m_arr[y][x]++;
					if (m_arr[y][x] >= LAM) flag1 = 1;
				}
				else
				{
					m_arr[(N + y + ay) % N][(N + x + ax) % N] += 2; phl_on_m_arr[(N + y + ay) % N][(N + x + ax) % N]++;
					if (m_arr[(N + y + ay) % N][(N + x + ax) % N] >= LAM) flag2 = 1;
				}
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////
void inits(void)         // Initialization of the system
{
	int j, m, k;
	seed(SD);

	grlength = 0;         // Caculate the length of the characteristic sequence of GR
	for (j = 0; grseq[j] != 0; j++)
		grlength++;

	nrlength = 0;         // Caculate the length of the characteristic sequence of NR
	for (j = 0; nrseq[j] != 0; j++)
		nrlength++;

	nprlength = 0;         // Caculate the length of the characteristic sequence of NPR
	for (j = 0; nprseq[j] != 0; j++)
		nprlength++;

	contrlength1 = 0;
	for (j = 0; contrseq1[j] != 0; j++)
		contrlength1++;

	contrlength2 = 0;
	for (j = 0; contrseq2[j] != 0; j++)
		contrlength2++;

	inoculength1 = 0;
	for (j = 0; inocuseq1[j] != 0; j++)
		inoculength1++;

	inoculength2 = 0;
	for (j = 0; inocuseq2[j] != 0; j++)
		inoculength2++;

	inoculength3 = 0;
	for (j = 0; inocuseq3[j] != 0; j++)
		inoculength3++;

	inoculength4 = 0;
	for (j = 0; inocuseq4[j] != 0; j++)
		inoculength4++;

	inoculength5 = 0;
	for (j = 0; inocuseq5[j] != 0; j++)
		inoculength5++;

	for (m = 0; m < 2; m++)     // Initiate the linked list array of nucleotides and RNA
	{
		for (y = 0; y < N; y++)
		{
			for (x = 0; x < N; x++)
			{
				p1 = (struct rna*)malloc(LEN);
				if (!p1) { printf("\tinit1--memeout\n"); exit(0); }
				room_head[m][y][x] = p1;
				p1->next = room_head[m][y][x];
			}
		}
	}

	for (k = 0; k < TNPPB; k++)  // Initial distribution of nucleotide-precursor's precursors
	{
		x = randl(N);
		y = randl(N);
		npp_arr[0][y][x]++;
	}

	for (k = 0; k < TFB; k++)   // Initial distribution of fatty acids
	{
		x = randl(N);
		y = randl(N);
		fa_arr[0][y][x]++;
	}

	for (k = 0; k < TGPB; k++)  // Initial distribution of glycerophosphate precursors
	{
		x = randl(N);
		y = randl(N);
		gp_arr[0][y][x]++;
	}
}

/////////////////////////////////////////////////////////////////////////
void inocu_cell(void)   // Inoculation of empty protocells
{
	int k;
	for (k = 0; k < INOCU_E_CELL_NUM; k++)
	{
		flag1 = 1;
		while (flag1)
		{
			x = randl(N);
			y = randl(N);
			if (m_arr[y][x] != 0)continue;
			m_arr[y][x] += LAM * 1.5;
			flag1 = 0;
		}
	}
}

void inoculate1(void)   // Inoculation of RNA
{
	int k, k1, k2;

	k = INOCU_CELL_NUM;
	for (y = 0; y < N; y++)   // Considering each room in the grid
	{
		if (k == 0)break;
		for (x = 0; x < N; x++)
		{
			if (k == 0)break;
			if (m_arr[y][x] > 0)  // Membrane existing (the room is occupied by a protocell)
			{
				k--;
				for (k1 = 0; k1 < INOCU_SEQ_NUM; k1++)
				{
					p2 = (struct rna*)malloc(LEN);
					if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
					for (k2 = 0; k2 < inoculength1; k2++) p2->information[0][k2] = inocuseq1[k2];
					p2->information[0][k2] = 0;
					p2->information[1][0] = 0;
					p2->length1 = inoculength1;
					p2->length2 = 0;
					p2->next = room_head[0][y][x]->next;
					if (p2->next != room_head[0][y][x])(p2->next)->prior = p2;
					room_head[0][y][x]->next = p2;
					p2->prior = room_head[0][y][x];
				}
			}
		}
	}
	if (k != 0)printf("no sufficient c for inoculating rna");

	k = INOCU_CELL_NUM;
	for (y = 0; y < N; y++)   // Considering each room in the grid
	{
		if (k == 0)break;
		for (x = 0; x < N; x++)
		{
			if (k == 0)break;
			if (m_arr[y][x] > 0)  // Membrane existing (the room is occupied by a protocell)
			{
				flag1 = 1;
				for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
				{
					if (findseq(inocuseq1, inoculength1, p) == 0) { flag1 = 0; break; }
				}
				if (flag1 != 0)
				{
					k--;
					for (k1 = 0; k1 < INOCU_SEQ_NUM; k1++)
					{
						p2 = (struct rna*)malloc(LEN);
						if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
						for (k2 = 0; k2 < inoculength3; k2++) p2->information[0][k2] = inocuseq3[k2];
						p2->information[0][k2] = 0;
						p2->information[1][0] = 0;
						p2->length1 = inoculength3;
						p2->length2 = 0;
						p2->next = room_head[0][y][x]->next;
						if (p2->next != room_head[0][y][x])(p2->next)->prior = p2;
						room_head[0][y][x]->next = p2;
						p2->prior = room_head[0][y][x];
					}
				}
			}
		}
	}
	if (k != 0) { printf("no sufficient c for inoculating rna");}
}

////////////////////////////////////////////////////////////////////////////////////////
void unit_action(void)      // Action (movement and events) of units (molecules and protocells) in the system
{
	int a, b, d, j, randnt, nr_turn, npr_turn, gr_turn, left, right, up, down, tmr_xy;
	double f, f1, f2, rtdaddlig, rtdaddphili, rntsyn, rnpsyn, rgsyn;

	//=================================== Movement of molecules
	//----------------- Fatty acids moving (including joining into membrane)
	for (y = 0; y < N; y++)   // Considering each room in the grid
	{
		for (x = 0; x < N; x++)
		{
			if (m_arr[y][x] > 0)  // Membrane existing (the room is occupied by a protocell)
			{					//  Possibly joining into the membrane
				fa_bef = fa_arr[0][y][x];
				fa_arr[0][y][x] = 0;
				for (k = 0; k < fa_bef; k++)
				{
					if (randd() < PMV && randd() < PFJM) m_arr[y][x]++;
					else fa_arr[1][y][x]++;
				}
			}
			else   // No membrane
			{		// Possibly moving, possibly joining into the membrane of a protocell at an adjacent room
				fa_bef = fa_arr[0][y][x];
				fa_arr[0][y][x] = 0;
				for (k = 0; k < fa_bef; k++)
				{
					if (randd() < PMV)
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:   // To left
							if (m_arr[y][(N + x - 1) % N] == 0) fa_arr[1][y][(N + x - 1) % N]++; // Moving
							else
							{
								if (randd() < PFJM) m_arr[y][(N + x - 1) % N]++; // Joining
								else fa_arr[1][y][x]++;
							}
							break;
						case 1:  // To right
							if (m_arr[y][(x + 1) % N] == 0) fa_arr[1][y][(x + 1) % N]++;
							else
							{
								if (randd() < PFJM) m_arr[y][(x + 1) % N]++;
								else fa_arr[1][y][x]++;
							}
							break;
						case 2:  // To up
							if (m_arr[(N + y - 1) % N][x] == 0) fa_arr[1][(N + y - 1) % N][x]++;
							else
							{
								if (randd() < PFJM) m_arr[(N + y - 1) % N][x]++;
								else fa_arr[1][y][x]++;
							}
							break;
						case 3:   // To down
							if (m_arr[(y + 1) % N][x] == 0) fa_arr[1][(y + 1) % N][x]++;
							else
							{
								if (randd() < PFJM) m_arr[(y + 1) % N][x]++;
								else fa_arr[1][y][x]++;
							}
							break;
						default: printf("fa moving error");
						}
					}
					else fa_arr[1][y][x]++;
				}
			}
		}
	}

	//--------------------------- Membrane components leaving membrane
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			if (m_arr[y][x] == 0)continue;  // No membrane
			tmr_xy = 0;
			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
				tmr_xy += p->length1 + p->length2; // Caculating total materials of nucleotides and RNA in the protocell
			// To consider the effect of osomitic pressure 
			m_fcd_bef = m_arr[y][x] - 2 * phl_on_m_arr[y][x];    // Fatty acids leaving
			for (k = 0; k < m_fcd_bef; k++)
			{
				if (randd() < PFLM / ((1 + tmr_xy / pow(m_arr[y][x] / 2.0, 1.5)) * (1 + FPL * 2.0 * phl_on_m_arr[y][x] / m_arr[y][x])))
				{                                          // Considering the influence of phospholipids on probability of fatty acids' leaving
					m_arr[y][x]--;  // Leaving
					if (randd() < 0.5) fa_arr[1][y][x]++;  // Going into the protocell
					else                              // Going out of the protocell 
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:
							if (m_arr[y][(N + x - 1) % N] == 0) fa_arr[1][y][(N + x - 1) % N]++;
							else m_arr[y][x]++;  // Blocked by adjacent membrane, cannot leave
							break;
						case 1:
							if (m_arr[y][(x + 1) % N] == 0) fa_arr[1][y][(x + 1) % N]++;
							else m_arr[y][x]++;
							break;
						case 2:
							if (m_arr[(N + y - 1) % N][x] == 0) fa_arr[1][(N + y - 1) % N][x]++;
							else m_arr[y][x]++;
							break;
						case 3:
							if (m_arr[(y + 1) % N][x] == 0) fa_arr[1][(y + 1) % N][x]++;
							else m_arr[y][x]++;
							break;
						default: printf("fatty acid leaving error");
						}
					}
				}
			}
			m_phl_bef = phl_on_m_arr[y][x];    // Phospholipids leaving
			for (k = 0; k < m_phl_bef; k++)
			{
				if (randd() < PPLM / ((1 + tmr_xy / pow(m_arr[y][x] / 2.0, 1.5)) * (1 + FPL * 2.0 * phl_on_m_arr[y][x] / m_arr[y][x])))
				{                                           // Considering the influence of phospholipids on probability of phospholipids' leaving
					m_arr[y][x] -= 2;  // Leaving
					phl_on_m_arr[y][x]--;
					if (randd() < 0.5) phl_arr[1][y][x]++;  // Going into the protocell
					else                              // Going out of the protocell 
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:
							if (m_arr[y][(N + x - 1) % N] == 0) phl_arr[1][y][(N + x - 1) % N]++;
							else { m_arr[y][x] += 2; phl_on_m_arr[y][x]++; } // Blocked by adjacent membrane, cannot leave
							break;
						case 1:
							if (m_arr[y][(x + 1) % N] == 0) phl_arr[1][y][(x + 1) % N]++;
							else { m_arr[y][x] += 2; phl_on_m_arr[y][x]++; }
							break;
						case 2:
							if (m_arr[(N + y - 1) % N][x] == 0) phl_arr[1][(N + y - 1) % N][x]++;
							else { m_arr[y][x] += 2; phl_on_m_arr[y][x]++; }
							break;
						case 3:
							if (m_arr[(y + 1) % N][x] == 0) phl_arr[1][(y + 1) % N][x]++;
							else { m_arr[y][x] += 2; phl_on_m_arr[y][x]++; }
							break;
						default: printf("phospholipid leaving error");
						}
					}
				}
			}
		}
	}

	//----------------- Phospholipids moving (including joining into membrane)
	for (y = 0; y < N; y++)   // Considering each room in the grid
	{
		for (x = 0; x < N; x++)
		{
			if (m_arr[y][x] > 0)  // Membrane existing (the room is occupied by a protocell)
			{					//  Possibly joining into the membrane
				phl_bef = phl_arr[0][y][x];
				phl_arr[0][y][x] = 0;
				for (k = 0; k < phl_bef; k++)
				{
					if (randd() < PMV && randd() < PPJM) { m_arr[y][x] += 2; phl_on_m_arr[y][x]++; }
					else phl_arr[1][y][x]++;
				}
			}
			else   // No membrane
			{		// Possibly moving, possibly joining into the membrane of a protocell at an adjacent room
				phl_bef = phl_arr[0][y][x];
				phl_arr[0][y][x] = 0;
				for (k = 0; k < phl_bef; k++)
				{
					if (randd() < PMV)
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:   // To left
							if (m_arr[y][(N + x - 1) % N] == 0) phl_arr[1][y][(N + x - 1) % N]++; // Moving
							else
							{
								if (randd() < PPJM) { m_arr[y][(N + x - 1) % N] += 2; phl_on_m_arr[y][(N + x - 1) % N]++; } // Joining
								else phl_arr[1][y][x]++;
							}
							break;
						case 1:  // To right
							if (m_arr[y][(x + 1) % N] == 0) phl_arr[1][y][(x + 1) % N]++;
							else
							{
								if (randd() < PPJM) { m_arr[y][(x + 1) % N] += 2; phl_on_m_arr[y][(x + 1) % N]++; }
								else phl_arr[1][y][x]++;
							}
							break;
						case 2:  // To up
							if (m_arr[(N + y - 1) % N][x] == 0) phl_arr[1][(N + y - 1) % N][x]++;
							else
							{
								if (randd() < PPJM) { m_arr[(N + y - 1) % N][x] += 2; phl_on_m_arr[(N + y - 1) % N][x]++; }
								else phl_arr[1][y][x]++;
							}
							break;
						case 3:   // To down
							if (m_arr[(y + 1) % N][x] == 0) phl_arr[1][(y + 1) % N][x]++;
							else
							{
								if (randd() < PPJM) { m_arr[(y + 1) % N][x] += 2; phl_on_m_arr[(y + 1) % N][x]++; }
								else phl_arr[1][y][x]++;
							}
							break;
						default: printf("pl moving error");
						}
					}
					else phl_arr[1][y][x]++;
				}
			}
		}
	}

	//-------------------Glycerophosphates moving	
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			if (m_arr[y][x] > 0) // Membrane existing
			{
				g_bef = g_arr[0][y][x];
				g_arr[0][y][x] = 0;
				for (k = 0; k < g_bef; k++)
				{
					if (randd() < PMV && randd() < PPF && (m_arr[y][x] - 2 * phl_on_m_arr[y][x]) >= 2) phl_on_m_arr[y][x]++; // Phospholipid forming
					else g_arr[1][y][x]++;  // No moving and no forming
				}
			}
			else   // No membrane
			{		// Possibly moving
				g_bef = g_arr[0][y][x];
				g_arr[0][y][x] = 0;
				for (k = 0; k < g_bef; k++)
				{
					if (randd() < PMV)
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:   // To left
							if (m_arr[y][(N + x - 1) % N] == 0) g_arr[1][y][(N + x - 1) % N]++; // Moving
							else if (randd() < PPF && (m_arr[y][(N + x - 1) % N] - 2 * phl_on_m_arr[y][(N + x - 1) % N]) >= 2) phl_on_m_arr[y][(N + x - 1) % N]++; // Phospholipid forming
							else g_arr[1][y][x]++;  // No moving and no forming
							break;
						case 1:  // To right
							if (m_arr[y][(x + 1) % N] == 0) g_arr[1][y][(x + 1) % N]++;
							else if (randd() < PPF && (m_arr[y][(x + 1) % N] - 2 * phl_on_m_arr[y][(x + 1) % N]) >= 2) phl_on_m_arr[y][(x + 1) % N]++; //** Phospholipid forming
							else g_arr[1][y][x]++;
							break;
						case 2:  // To up
							if (m_arr[(N + y - 1) % N][x] == 0) g_arr[1][(N + y - 1) % N][x]++;
							else if (randd() < PPF && (m_arr[(N + y - 1) % N][x] - 2 * phl_on_m_arr[(N + y - 1) % N][x]) >= 2) phl_on_m_arr[(N + y - 1) % N][x]++; //** Phospholipid forming								
							else g_arr[1][y][x]++;
							break;
						case 3:   // To down
							if (m_arr[(y + 1) % N][x] == 0) g_arr[1][(y + 1) % N][x]++;
							else if (randd() < PPF && (m_arr[(y + 1) % N][x] - 2 * phl_on_m_arr[(y + 1) % N][x]) >= 2) phl_on_m_arr[(y + 1) % N][x]++; //** Phospholipid forming
							else g_arr[1][y][x]++;
							break;
						default: printf("g moving error");
						}
					}
					else g_arr[1][y][x]++;
				}
			}
		}
	}

	//------------------------------ Nucleotides and RNA moving
	if (i < MECH_CHANGE_STEP) f2 = 0;
	else f2 = FPP;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			if (m_arr[y][x] > 0) // Membrane existing
			{
				for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
				{
					if (p->length1 + p->length2 > 1) fresh_rna(0);  // Non-nt cannot permeate
					else
					{
						if (randd() < PMV && randd() < PNP / (1 + f2 * 2.0 * phl_on_m_arr[y][x] / m_arr[y][x]))   //Nucleotides permeating               
						{
							randcase = randl(4);   // Four possible directions
							switch (randcase)
							{
							case 0:
								if (m_arr[y][(N + x - 1) % N] == 0)
								{
									p1 = p->prior;
									p2 = p->next;
									p3 = room_head[1][y][(N + x - 1) % N]->next;
									room_head[1][y][(N + x - 1) % N]->next = p;
									p->next = p3;
									p->prior = room_head[1][y][(N + x - 1) % N];
									if (p3 != room_head[1][y][(N + x - 1) % N])p3->prior = p;
									p1->next = p2;
									if (p2 != room_head[0][y][x])p2->prior = p1;
									p = p1;
								}
								else fresh_rna(0);
								break;
							case 1:
								if (m_arr[y][(x + 1) % N] == 0)
								{
									p1 = p->prior;
									p2 = p->next;
									p3 = room_head[1][y][(x + 1) % N]->next;
									room_head[1][y][(x + 1) % N]->next = p;
									p->next = p3;
									p->prior = room_head[1][y][(x + 1) % N];
									if (p3 != room_head[1][y][(x + 1) % N])p3->prior = p;
									p1->next = p2;
									if (p2 != room_head[0][y][x])p2->prior = p1;
									p = p1;
								}
								else fresh_rna(0);
								break;
							case 2:
								if (m_arr[(N + y - 1) % N][x] == 0)
								{
									p1 = p->prior;
									p2 = p->next;
									p3 = room_head[1][(N + y - 1) % N][x]->next;
									room_head[1][(N + y - 1) % N][x]->next = p;
									p->next = p3;
									p->prior = room_head[1][(N + y - 1) % N][x];
									if (p3 != room_head[1][(N + y - 1) % N][x])p3->prior = p;
									p1->next = p2;
									if (p2 != room_head[0][y][x])p2->prior = p1;
									p = p1;
								}
								else fresh_rna(0);
								break;
							case 3:
								if (m_arr[(y + 1) % N][x] == 0)
								{
									p1 = p->prior;
									p2 = p->next;
									p3 = room_head[1][(y + 1) % N][x]->next;
									room_head[1][(y + 1) % N][x]->next = p;
									p->next = p3;
									p->prior = room_head[1][(y + 1) % N][x];
									if (p3 != room_head[1][(y + 1) % N][x])p3->prior = p;
									p1->next = p2;
									if (p2 != room_head[0][y][x])p2->prior = p1;
									p = p1;
								}
								else fresh_rna(0);
								break;
							default: printf("rna moving error");
							}
						}
						else fresh_rna(0);
					}
				}
			}
			else   // No membrane
			{
				for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
				{
					if (randd() * RMRW < PMV)
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:
							if (m_arr[y][(N + x - 1) % N] == 0)  // No membrane existing in this direction
							{
								p1 = p->prior;
								p2 = p->next;
								p3 = room_head[1][y][(N + x - 1) % N]->next;
								room_head[1][y][(N + x - 1) % N]->next = p;
								p->next = p3;
								p->prior = room_head[1][y][(N + x - 1) % N];
								if (p3 != room_head[1][y][(N + x - 1) % N])p3->prior = p;
								p1->next = p2;
								if (p2 != room_head[0][y][x])p2->prior = p1;
								p = p1;
							}
							else   // Membrane existing in this direction
							{
								if (p->length1 + p->length2 > 1) fresh_rna(0);  // Non_nt cannot permeate
								else        //Nucleotides permeating
								{
									tmr_xy = 0;
									for (p4 = room_head[0][y][(N + x - 1) % N]->next; p4 != room_head[0][y][(N + x - 1) % N]; p4 = p4->next)
									{
										if (p4->length1 + p4->length2 > 1)
										{
											tmr_xy += p4->length1 + p4->length2; // Caculating total materials of non-nt RNA (thus impermeable) in the target protocell 
										}
									}
									// To consider the effect of Donnan's equilibrium 	
									if (randd() < PNP * (m_arr[y][(N + x - 1) % N] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[y][(N + x - 1) % N] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[y][(N + x - 1) % N] / m_arr[y][(N + x - 1) % N])))
									{
										p1 = p->prior;
										p2 = p->next;
										p3 = room_head[1][y][(N + x - 1) % N]->next;
										room_head[1][y][(N + x - 1) % N]->next = p;
										p->next = p3;
										p->prior = room_head[1][y][(N + x - 1) % N];
										if (p3 != room_head[1][y][(N + x - 1) % N])p3->prior = p;
										p1->next = p2;
										if (p2 != room_head[0][y][x])p2->prior = p1;
										p = p1;
									}
									else fresh_rna(0);
								}
							}
							break;
						case 1:
							if (m_arr[y][(x + 1) % N] == 0)  // No membrane existing in this direction
							{
								p1 = p->prior;
								p2 = p->next;
								p3 = room_head[1][y][(x + 1) % N]->next;
								room_head[1][y][(x + 1) % N]->next = p;
								p->next = p3;
								p->prior = room_head[1][y][(x + 1) % N];
								if (p3 != room_head[1][y][(x + 1) % N])p3->prior = p;
								p1->next = p2;
								if (p2 != room_head[0][y][x])p2->prior = p1;
								p = p1;
							}
							else   // Membrane existing in this direction
							{
								if (p->length1 + p->length2 > 1) fresh_rna(0);  // Non_nt cannot permeate
								else               //Nucleotides permeating
								{
									tmr_xy = 0;
									for (p4 = room_head[0][y][(x + 1) % N]->next; p4 != room_head[0][y][(x + 1) % N]; p4 = p4->next)
									{
										if (p4->length1 + p4->length2 > 1)
										{
											tmr_xy += p4->length1 + p4->length2; // Caculating total materials of nucleotides and RNA in the target protocell 
										}
									}
									// To consider the effect of Donnan's equilibrium 	
									if (randd() < PNP * (m_arr[y][(x + 1) % N] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[y][(x + 1) % N] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[y][(x + 1) % N] / m_arr[y][(x + 1) % N])))
									{
										p1 = p->prior;
										p2 = p->next;
										p3 = room_head[1][y][(x + 1) % N]->next;
										room_head[1][y][(x + 1) % N]->next = p;
										p->next = p3;
										p->prior = room_head[1][y][(x + 1) % N];
										if (p3 != room_head[1][y][(x + 1) % N])p3->prior = p;
										p1->next = p2;
										if (p2 != room_head[0][y][x])p2->prior = p1;
										p = p1;
									}
									else fresh_rna(0);
								}
							}
							break;
						case 2:
							if (m_arr[(N + y - 1) % N][x] == 0)
							{
								p1 = p->prior;
								p2 = p->next;
								p3 = room_head[1][(N + y - 1) % N][x]->next;
								room_head[1][(N + y - 1) % N][x]->next = p;
								p->next = p3;
								p->prior = room_head[1][(N + y - 1) % N][x];
								if (p3 != room_head[1][(N + y - 1) % N][x])p3->prior = p;
								p1->next = p2;
								if (p2 != room_head[0][y][x])p2->prior = p1;
								p = p1;
							}
							else
							{
								if (p->length1 + p->length2 > 1) fresh_rna(0);
								else
								{
									tmr_xy = 0;
									for (p4 = room_head[0][(N + y - 1) % N][x]->next; p4 != room_head[0][(N + y - 1) % N][x]; p4 = p4->next)
									{
										if (p4->length1 + p4->length2 > 1)
										{
											tmr_xy += p4->length1 + p4->length2;
										}
									}
									if (randd() < PNP * (m_arr[(N + y - 1) % N][x] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[(N + y - 1) % N][x] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[(N + y - 1) % N][x] / m_arr[(N + y - 1) % N][x])))
									{
										p1 = p->prior;
										p2 = p->next;
										p3 = room_head[1][(N + y - 1) % N][x]->next;
										room_head[1][(N + y - 1) % N][x]->next = p;
										p->next = p3;
										p->prior = room_head[1][(N + y - 1) % N][x];
										if (p3 != room_head[1][(N + y - 1) % N][x])p3->prior = p;
										p1->next = p2;
										if (p2 != room_head[0][y][x])p2->prior = p1;
										p = p1;
									}
									else fresh_rna(0);
								}
							}
							break;
						case 3:
							if (m_arr[(y + 1) % N][x] == 0)
							{
								p1 = p->prior;
								p2 = p->next;
								p3 = room_head[1][(y + 1) % N][x]->next;
								room_head[1][(y + 1) % N][x]->next = p;
								p->next = p3;
								p->prior = room_head[1][(y + 1) % N][x];
								if (p3 != room_head[1][(y + 1) % N][x])p3->prior = p;
								p1->next = p2;
								if (p2 != room_head[0][y][x])p2->prior = p1;
								p = p1;
							}
							else
							{
								if (p->length1 + p->length2 > 1) fresh_rna(0);
								else
								{
									tmr_xy = 0;
									for (p4 = room_head[0][(y + 1) % N][x]->next; p4 != room_head[0][(y + 1) % N][x]; p4 = p4->next)
									{
										if (p4->length1 + p4->length2 > 1)
										{
											tmr_xy += p4->length1 + p4->length2;
										}
									}

									if (randd() < PNP * (m_arr[(y + 1) % N][x] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[(y + 1) % N][x] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[(y + 1) % N][x] / m_arr[(y + 1) % N][x])))
									{
										p1 = p->prior;
										p2 = p->next;
										p3 = room_head[1][(y + 1) % N][x]->next;
										room_head[1][(y + 1) % N][x]->next = p;
										p->next = p3;
										p->prior = room_head[1][(y + 1) % N][x];
										if (p3 != room_head[1][(y + 1) % N][x])p3->prior = p;
										p1->next = p2;
										if (p2 != room_head[0][y][x])p2->prior = p1;
										p = p1;
									}
									else fresh_rna(0);
								}
							}
							break;
						default: printf("rna moving error");
						}
					}
					else fresh_rna(0);
				}
			}
		}
	}

	//-----------------Glycerophosphate precursors moving	
	if (i < MECH_CHANGE_STEP) f2 = 0;
	else f2 = FPPW;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			if (m_arr[y][x] > 0)  // Membrane existing
			{					// Possibly moving and permeating out
				gp_bef = gp_arr[0][y][x];
				gp_arr[0][y][x] = 0;
				for (k = 0; k < gp_bef; k++)
				{
					if (randd() < PMV && randd() < PGPP / (1 + f2 * 2.0 * phl_on_m_arr[y][x] / m_arr[y][x]))
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:
							if (m_arr[y][(N + x - 1) % N] == 0) gp_arr[1][y][(N + x - 1) % N]++;
							else gp_arr[1][y][x]++;
							break;
						case 1:
							if (m_arr[y][(x + 1) % N] == 0) gp_arr[1][y][(x + 1) % N]++;
							else gp_arr[1][y][x]++;
							break;
						case 2:
							if (m_arr[(N + y - 1) % N][x] == 0) gp_arr[1][(N + y - 1) % N][x]++;
							else gp_arr[1][y][x]++;
							break;
						case 3:
							if (m_arr[(y + 1) % N][x] == 0) gp_arr[1][(y + 1) % N][x]++;
							else gp_arr[1][y][x]++;
							break;
						default: printf("gp moving error");
						}
					}
					else gp_arr[1][y][x]++;
				}
			}
			else    // No membrane
			{		// Possibly moving, possibly permeating into a protocell at an adjacent room 
				gp_bef = gp_arr[0][y][x];
				gp_arr[0][y][x] = 0;
				for (k = 0; k < gp_bef; k++)
				{
					if (randd() < PMV)
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:
							if (m_arr[y][(N + x - 1) % N] == 0) gp_arr[1][y][(N + x - 1) % N]++;   // No membrane existing in this direction
							else
							{
								if (randd() < PGPP * (m_arr[y][(N + x - 1) % N] / LAM) / (1 + f2 * 2.0 * phl_on_m_arr[y][(N + x - 1) % N] / m_arr[y][(N + x - 1) % N]))
									gp_arr[1][y][(N + x - 1) % N]++;
								else gp_arr[1][y][x]++;
							}
							break;
						case 1:
							if (m_arr[y][(x + 1) % N] == 0) gp_arr[1][y][(x + 1) % N]++;
							else
							{
								if (randd() < PGPP * (m_arr[y][(x + 1) % N] / LAM) / (1 + f2 * 2.0 * phl_on_m_arr[y][(x + 1) % N] / m_arr[y][(x + 1) % N]))
									gp_arr[1][y][(x + 1) % N]++;
								else gp_arr[1][y][x]++;
							}
							break;
						case 2:
							if (m_arr[(N + y - 1) % N][x] == 0) gp_arr[1][(N + y - 1) % N][x]++;
							else
							{
								if (randd() < PGPP * (m_arr[(N + y - 1) % N][x] / LAM) / (1 + f2 * 2.0 * phl_on_m_arr[(N + y - 1) % N][x] / m_arr[(N + y - 1) % N][x]))
									gp_arr[1][(N + y - 1) % N][x]++;
								else gp_arr[1][y][x]++;
							}
							break;
						case 3:
							if (m_arr[(y + 1) % N][x] == 0) gp_arr[1][(y + 1) % N][x]++;
							else
							{
								if (randd() < PGPP * (m_arr[(y + 1) % N][x] / LAM) / (1 + f2 * 2.0 * phl_on_m_arr[(y + 1) % N][x] / m_arr[(y + 1) % N][x]))
									gp_arr[1][(y + 1) % N][x]++;
								else gp_arr[1][y][x]++;
							}
							break;
						default: printf("pl moving error");
						}
					}
					else gp_arr[1][y][x]++;
				}
			}
		}
	}

	//--------------------------------- Nucleotide precursors moving
	if (i < MECH_CHANGE_STEP) f2 = 0;
	else f2 = FPP;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			if (m_arr[y][x] > 0)  // Membrane existing
			{				   // Possibly moving and permeating out
				np_bef = np_arr[0][y][x];
				np_arr[0][y][x] = 0;
				for (k = 0; k < np_bef; k++)
				{
					if (randd() < PMV && randd() < PNPP / (1 + f2 * 2.0 * phl_on_m_arr[y][x] / m_arr[y][x]))
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:
							if (m_arr[y][(N + x - 1) % N] == 0) np_arr[1][y][(N + x - 1) % N]++;
							else np_arr[1][y][x]++;
							break;
						case 1:
							if (m_arr[y][(x + 1) % N] == 0) np_arr[1][y][(x + 1) % N]++;
							else np_arr[1][y][x]++;
							break;
						case 2:
							if (m_arr[(N + y - 1) % N][x] == 0) np_arr[1][(N + y - 1) % N][x]++;
							else np_arr[1][y][x]++;
							break;
						case 3:
							if (m_arr[(y + 1) % N][x] == 0) np_arr[1][(y + 1) % N][x]++;
							else np_arr[1][y][x]++;
							break;
						default: printf("np moving error");
						}
					}
					else np_arr[1][y][x]++;
				}
			}
			else    // no membrane
			{		// Possible moving, possibly permeating into a protocell at an adjacent room
				np_bef = np_arr[0][y][x];
				np_arr[0][y][x] = 0;
				for (k = 0; k < np_bef; k++)
				{
					if (randd() < PMV)
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:
							if (m_arr[y][(N + x - 1) % N] == 0) np_arr[1][y][(N + x - 1) % N]++; // Moving
							else
							{
								tmr_xy = 0;
								for (p = room_head[1][y][(N + x - 1) % N]->next; p != room_head[1][y][(N + x - 1) % N]; p = p->next)
								{
									if (p->length1 + p->length2 > 1)
									{
										tmr_xy += p->length1 + p->length2; // Caculating total materials of non-nt RNA (thus impermeable) in the target protocell 
									}
								}
								// To consider the effect of Donnan's equilibrium 	
								if (randd() < PNPP * (m_arr[y][(N + x - 1) % N] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[y][(N + x - 1) % N] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[y][(N + x - 1) % N] / m_arr[y][(N + x - 1) % N])))
									np_arr[1][y][(N + x - 1) % N]++;
								else np_arr[1][y][x]++;
							}
							break;
						case 1:
							if (m_arr[y][(x + 1) % N] == 0) np_arr[1][y][(x + 1) % N]++;
							else
							{
								tmr_xy = 0;
								for (p = room_head[1][y][(x + 1) % N]->next; p != room_head[1][y][(x + 1) % N]; p = p->next)
								{
									if (p->length1 + p->length2 > 1)
									{
										tmr_xy += p->length1 + p->length2;
									}
								}
								if (randd() < PNPP * (m_arr[y][(x + 1) % N] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[y][(x + 1) % N] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[y][(x + 1) % N] / m_arr[y][(x + 1) % N])))
									np_arr[1][y][(x + 1) % N]++;
								else np_arr[1][y][x]++;
							}
							break;
						case 2:
							if (m_arr[(N + y - 1) % N][x] == 0) np_arr[1][(N + y - 1) % N][x]++;
							else
							{
								tmr_xy = 0;
								for (p = room_head[1][(N + y - 1) % N][x]->next; p != room_head[1][(N + y - 1) % N][x]; p = p->next)
								{
									if (p->length1 + p->length2 > 1)
									{
										tmr_xy += p->length1 + p->length2;
									}
								}
								if (randd() < PNPP * (m_arr[(N + y - 1) % N][x] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[(N + y - 1) % N][x] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[(N + y - 1) % N][x] / m_arr[(N + y - 1) % N][x])))
									np_arr[1][(N + y - 1) % N][x]++;
								else np_arr[1][y][x]++;
							}
							break;
						case 3:
							if (m_arr[(y + 1) % N][x] == 0) np_arr[1][(y + 1) % N][x]++;
							else
							{
								tmr_xy = 0;
								for (p = room_head[1][(y + 1) % N][x]->next; p != room_head[1][(y + 1) % N][x]; p = p->next)
								{
									if (p->length1 + p->length2 > 1)
									{
										tmr_xy += p->length1 + p->length2;
									}
								}
								if (randd() < PNPP * (m_arr[(y + 1) % N][x] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[(y + 1) % N][x] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[(y + 1) % N][x] / m_arr[(y + 1) % N][x])))
									np_arr[1][(y + 1) % N][x]++;
								else np_arr[1][y][x]++;
							}
							break;
						default: printf("np moving error");
						}
					}
					else np_arr[1][y][x]++;
				}
			}
		}
	}

	//--------------------------------- Nucleotide precursors' precursors moving
	if (i < MECH_CHANGE_STEP) f2 = 0;
	else f2 = FPPW;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			if (m_arr[y][x] > 0)  // Membrane existing
			{				   // Possibly moving and permeating out
				npp_bef = npp_arr[0][y][x];
				npp_arr[0][y][x] = 0;
				for (k = 0; k < npp_bef; k++)
				{
					if (randd() < PMV && randd() < PNPPP / (1 + f2 * 2.0 * phl_on_m_arr[y][x] / m_arr[y][x]))
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:
							if (m_arr[y][(N + x - 1) % N] == 0) npp_arr[1][y][(N + x - 1) % N]++;
							else npp_arr[1][y][x]++;
							break;
						case 1:
							if (m_arr[y][(x + 1) % N] == 0) npp_arr[1][y][(x + 1) % N]++;
							else npp_arr[1][y][x]++;
							break;
						case 2:
							if (m_arr[(N + y - 1) % N][x] == 0) npp_arr[1][(N + y - 1) % N][x]++;
							else npp_arr[1][y][x]++;
							break;
						case 3:
							if (m_arr[(y + 1) % N][x] == 0) npp_arr[1][(y + 1) % N][x]++;
							else npp_arr[1][y][x]++;
							break;
						default: printf("npp moving error");
						}
					}
					else npp_arr[1][y][x]++;
				}
			}
			else    // No membrane
			{		// Possible moving, possibly permeating into a protocell at an adjacent room
				npp_bef = npp_arr[0][y][x];
				npp_arr[0][y][x] = 0;
				for (k = 0; k < npp_bef; k++)
				{
					if (randd() < PMV)
					{
						randcase = randl(4);   // Four possible directions
						switch (randcase)
						{
						case 0:
							if (m_arr[y][(N + x - 1) % N] == 0) npp_arr[1][y][(N + x - 1) % N]++; // Moving
							else
							{
								tmr_xy = 0;
								for (p = room_head[1][y][(N + x - 1) % N]->next; p != room_head[1][y][(N + x - 1) % N]; p = p->next)
								{
									if (p->length1 + p->length2 > 1)
									{
										tmr_xy += p->length1 + p->length2; // Caculating total materials of non-nt RNA (thus impermeable) in the target protocell 
									}
								}
								// To consider the effect of Donnan's equilibrium 	
								if (randd() < PNPPP * (m_arr[y][(N + x - 1) % N] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[y][(N + x - 1) % N] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[y][(N + x - 1) % N] / m_arr[y][(N + x - 1) % N])))
									npp_arr[1][y][(N + x - 1) % N]++;
								else npp_arr[1][y][x]++;
							}
							break;
						case 1:
							if (m_arr[y][(x + 1) % N] == 0) npp_arr[1][y][(x + 1) % N]++;
							else
							{
								tmr_xy = 0;
								for (p = room_head[1][y][(x + 1) % N]->next; p != room_head[1][y][(x + 1) % N]; p = p->next)
								{
									if (p->length1 + p->length2 > 1)
									{
										tmr_xy += p->length1 + p->length2;
									}
								}
								if (randd() < PNPPP * (m_arr[y][(x + 1) % N] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[y][(x + 1) % N] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[y][(x + 1) % N] / m_arr[y][(x + 1) % N])))
									npp_arr[1][y][(x + 1) % N]++;
								else npp_arr[1][y][x]++;
							}
							break;
						case 2:
							if (m_arr[(N + y - 1) % N][x] == 0) npp_arr[1][(N + y - 1) % N][x]++;
							else
							{
								tmr_xy = 0;
								for (p = room_head[1][(N + y - 1) % N][x]->next; p != room_head[1][(N + y - 1) % N][x]; p = p->next)
								{
									if (p->length1 + p->length2 > 1)
									{
										tmr_xy += p->length1 + p->length2;
									}
								}
								if (randd() < PNPPP * (m_arr[(N + y - 1) % N][x] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[(N + y - 1) % N][x] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[(N + y - 1) % N][x] / m_arr[(N + y - 1) % N][x])))
									npp_arr[1][(N + y - 1) % N][x]++;
								else npp_arr[1][y][x]++;
							}
							break;
						case 3:
							if (m_arr[(y + 1) % N][x] == 0) npp_arr[1][(y + 1) % N][x]++;
							else
							{
								tmr_xy = 0;
								for (p = room_head[1][(y + 1) % N][x]->next; p != room_head[1][(y + 1) % N][x]; p = p->next)
								{
									if (p->length1 + p->length2 > 1)
									{
										tmr_xy += p->length1 + p->length2;
									}
								}
								if (randd() < PNPPP * (m_arr[(y + 1) % N][x] / LAM) / ((1.0 + FDE * tmr_xy / pow(m_arr[(y + 1) % N][x] / 2.0, 1.5)) * (1 + f2 * 2.0 * phl_on_m_arr[(y + 1) % N][x] / m_arr[(y + 1) % N][x])))
									npp_arr[1][(y + 1) % N][x]++;
								else npp_arr[1][y][x]++;
							}
							break;
						default: printf("npp moving error");
						}
					}
					else npp_arr[1][y][x]++;
				}
			}
		}
	}
	//====================================== End of "Movement of molecules"

	//============================================= Events of molecules
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			//----------------------------Fatty acids' events
			fa_arr[0][y][x] = fa_arr[1][y][x];   // No events
			fa_arr[1][y][x] = 0;

			//--------------------------- phospholipids' events
			phl_bef = phl_arr[1][y][x];  // phospholipids out of membrane decaying
			phl_arr[1][y][x] = 0;
			for (k = 0; k < phl_bef; k++)
			{
				if (randd() < PPD) { g_arr[0][y][x]++; fa_arr[0][y][x] += 2; }
				else phl_arr[0][y][x]++;
			}
			if (m_arr[y][x] > 0)    // phospholipids on membrane decaying
			{
				m_phl_bef = phl_on_m_arr[y][x];
				for (k = 0; k < m_phl_bef; k++)
				{
					if (randd() < PPDM)
					{
						phl_on_m_arr[y][x]--;
						if (randd() < 0.5) g_arr[0][y][x]++;  // glycerophosphate going into the protocell
						else							// glycerophosphate going out of the protocell
						{
							randcase = randl(4);
							switch (randcase)
							{
							case 0:
								if (m_arr[y][(N + x - 1) % N] == 0) g_arr[0][y][(N + x - 1) % N]++;
								else phl_on_m_arr[y][x]++;   // There is a cell in this direction, thus no decay.
								break;
							case 1:
								if (m_arr[y][(x + 1) % N] == 0) g_arr[0][y][(x + 1) % N]++;
								else phl_on_m_arr[y][x]++;
								break;
							case 2:
								if (m_arr[(N + y - 1) % N][x] == 0) g_arr[0][(N + y - 1) % N][x]++;
								else phl_on_m_arr[y][x]++;
								break;
							case 3:
								if (m_arr[(y + 1) % N][x] == 0) g_arr[0][(y + 1) % N][x]++;
								else phl_on_m_arr[y][x]++;
								break;
							default: printf("phospholipids on membrane decaying error");
							}
						}
					}
				}
			}

			//---------------------- Glycerophosphates' events
			g_bef = g_arr[1][y][x];  // Glycerophosphate decaying
			g_arr[1][y][x] = 0;
			for (k = 0; k < g_bef; k++)
			{
				if (randd() < PGD) gp_arr[0][y][x]++;
				else g_arr[0][y][x]++;
			}

			//------------------------------ Nucleotide precurors' events
			np_bef = np_arr[1][y][x];  // Nucleotide precuror decaying 
			np_arr[1][y][x] = 0;
			for (k = 0; k < np_bef; k++)
			{
				if (randd() < PNPD) npp_arr[0][y][x]++;
				else np_arr[1][y][x]++;      // If not decaying, left for consideration of forming nt below
			}

			//------------------------------ Nucleotides and RNA's events
			rna_shuffle(1); // Random order-arrangement of nodes in the linked list for nucleotides and RNA in the grid room 
			for (p = room_head[1][y][x]->next; p != room_head[1][y][x]; p = p->next)
			{
				// Random chain ligating 
				for (p3 = p->next; p3 != p; p3 = p3->next)
				{
					if (p3 == room_head[1][y][x]) { p3 = room_head[1][y][x]->next; if (p3 == p)break; }
					if (p3->length2 == 0)
					{
						if (randd() < PRL / (p->length1 * p3->length1)) // Random ligating should be more difficult for longer chains
						{
							if (p->length1 + p3->length1 > MAX_RNA_LENGTH - 1) { over_max_len++; continue; }
							for (a = 0; a < p3->length1; a++) p->information[0][a + p->length1] = p3->information[0][a];
							p->information[0][p->length1 + p3->length1] = 0;
							p->length1 = p->length1 + p3->length1;
							(p3->prior)->next = p3->next;
							if (p3->next != room_head[1][y][x])(p3->next)->prior = p3->prior;
							free(p3);
							break;
						}
					}
				}

				// Decaying and degradating 
				if (p->length1 == 1)  // Nucleotide decaying
				{
					if (m_arr[y][x] == 0) f = PND * FDO;
					else f = PND;
					if (p->length2 == 0 && randd() < f)
					{
						np_arr[0][y][x]++;
						(p->prior)->next = p->next;
						if (p->next != room_head[1][y][x])(p->next)->prior = p->prior;
						p3 = p;
						p = p->prior;
						free(p3); continue;//break;
					}
				}
				else               // RNA end-decaying and degradating 
				{
					if (m_arr[y][x] == 0) f = PNDE * FDO;
					else f = PNDE;
					if (p->length1 > p->length2 && randd() < f) // Nucleotide residue decaying at the end of an RNA
					{
						np_arr[0][y][x]++;
						p->information[0][p->length1 - 1] = 0;
						p->length1--;
					}

					if (p->length1 != 1)    // RNA degradating
					{
						if (m_arr[y][x] == 0) f = PBB * FDO;
						else f = PBB;
						for (j = p->length1; j > 1; j--)
						{
							if (j <= p->length2 && p->nick[j - 1] == 0) f1 = f * sqrt(f);  // Breaking of double chain should be more difficult
							else f1 = f;

							if (randd() < f1)
							{
								p3 = (struct rna*)malloc(LEN);
								if (!p3) { printf("\t%ddeg--memeout\n", i); exit(0); }

								for (b = 0; b < p->length1 - j + 1; b++) p3->information[0][b] = p->information[0][b + j - 1];
								p3->information[0][p->length1 - j + 1] = 0;
								p->information[0][j - 1] = 0;
								p3->length1 = p->length1 - j + 1;
								p->length1 = j - 1;

								if (p->length2 > j - 1)
								{
									for (b = 0; b < p->length2 - j + 1; b++)	p3->information[1][b] = p->information[1][b + j - 1];
									p3->information[1][p->length2 - j + 1] = 0;
									p->information[1][j - 1] = 0;
									p3->length2 = p->length2 - j + 1;
									p->length2 = j - 1;
								}
								else
								{
									p3->information[1][0] = 0;
									p3->length2 = 0;
								}

								if (p3->length2 != 0)
								{
									p3->nick[0] = 1;
									for (a = 1; a < p3->length2; a++)
										p3->nick[a] = p->nick[j + a - 1];
								}

								p3->prior = room_head[1][y][x];
								p3->next = room_head[1][y][x]->next;
								if (p3->next != room_head[1][y][x])(p3->next)->prior = p3;
								room_head[1][y][x]->next = p3;
								break;
							}
						}
					}
				}

				//Template-directed attracting
				if (p->length1 > 3 && (p->length2 != 0 || randd() < 0.5)) // p->length1>3: If an RNA is shorter than 3 nt, it would not serve as a template
				{                                          // randd()<0.5: Single-chain RNA turning into template
					for (p3 = p->next; p3 != p; p3 = p3->next)   //Template-directed attraction of substrates
					{
						if (p3 == room_head[1][y][x]) { p3 = room_head[1][y][x]->next; if (p3 == p)break; }
						if (p3->length2 == 0)
						{
							if (p3->length1 <= p->length1 - p->length2)
							{
								for (flag1 = 0, b = 0; b < p3->length1; b++)
								{
									if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][p->length2 + b]) == 5)continue;
									else if (randd() < PFP)continue;
									else { flag1 = 1; break; }
								}
								if (flag1 == 0)
								{
									rtdaddphili = randd();
									if (rtdaddphili < PAT)
									{
										for (a = 0; a < p3->length1; a++)
											p->information[1][p->length2 + a] = p3->information[0][p3->length1 - 1 - a];
										p->information[1][p->length2 + p3->length1] = 0;

										p->nick[p->length2] = 1;
										for (a = 1; a < p3->length1; a++)
											p->nick[p->length2 + a] = 0;
										p->length2 = p->length2 + p3->length1;

										(p3->prior)->next = p3->next;
										if (p3->next != room_head[1][y][x])(p3->next)->prior = p3->prior;
										free(p3);
										break;
									}
								}
							}
						}
					}
				}

				//Template-directed ligating
				for (a = 1; a < p->length2; a++)
				{
					if (p->nick[a] == 0) continue;
					rtdaddlig = randd();
					if (rtdaddlig < PTL) p->nick[a] = 0;
				}

				// Double chain separating
				if (p->length2 != 0)
				{
					a = p->length2 - 1;
					while (p->nick[a] == 0) a--;

					if (randd() < pow(PSP, sqrt(1.0 * (p->length2 - a))))     // sqrt: considering single-chain folding factor    
					{
						p3 = (struct rna*)malloc(LEN);
						if (!p3) { printf("\t%dsep--memeout\n", i); exit(0); }

						for (b = 0; b < p->length2 - a; b++)
							p3->information[0][b] = p->information[1][p->length2 - 1 - b];
						p->information[1][a] = 0;

						p3->information[0][p->length2 - a] = 0;
						p3->information[1][0] = 0;

						p3->length1 = p->length2 - a;
						p->length2 = a;

						p3->length2 = 0;

						p3->prior = room_head[1][y][x];
						p3->next = room_head[1][y][x]->next;
						if (p3->next != room_head[1][y][x])(p3->next)->prior = p3;
						room_head[1][y][x]->next = p3;
					}
				}
			}

			for (p = room_head[1][y][x]->next; p != room_head[1][y][x]; p = p->next) fresh_rna(1);

			//---------------------------- Nucleotide-precuror's precursors' events
			npr_turn = 0;
			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
			{									// Finding NPR in the room
				if (p->length1 > 1.5 * nprlength)continue;   // Only an RNA no longer than 1.5 times of the characteristic sequence's length can act as the ribozyme
				flag1 = findseq(nprseq, nprlength, p);
				if (flag1 == 0 && p->length2 == 0) npr_turn += NPR_TURN;   // Caculating possible catalytic turns of the ribozyme in the room 
			}
			npp_bef = npp_arr[1][y][x];
			npp_arr[1][y][x] = 0;
			for (k = 0; k < npp_bef; k++)
			{
				rnpsyn = randd();
				flagnpsyn = 1;
				if (rnpsyn < PNPF) flagnpsyn = 0; // Uncatalyzed
				else if (npr_turn > 0)
				{
					npr_turn--;
					if (rnpsyn < PNPFR) flagnpsyn = 0; // Catalyzed
				}
				if (flagnpsyn == 0) np_arr[0][y][x]++; // Forming a nucleotide precursor 
				else npp_arr[0][y][x]++;
			}

			//------------------------------ Glycerophosphate precursors' events
			gr_turn = 0;
			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
			{									// Finding GR in the room
				if (p->length1 > 1.5 * grlength)continue;
				flag1 = findseq(grseq, grlength, p);
				if (flag1 == 0 && p->length2 == 0) gr_turn += GR_TURN;
			}
			gp_bef = gp_arr[1][y][x];
			gp_arr[1][y][x] = 0;
			for (k = 0; k < gp_bef; k++)
			{
				rgsyn = randd();
				flaggsyn = 1;
				if (rgsyn < PGF) flaggsyn = 0; // Uncatalyzed
				else if (gr_turn > 0)
				{
					gr_turn--;
					if (rgsyn < PGFR) flaggsyn = 0; // Catalyzed
				}
				if (flaggsyn == 0) g_arr[0][y][x]++; // Forming a glycerophosphate
				else gp_arr[0][y][x]++;
			}

			//--------------------------------  Nucleotide precursors' events
			nr_turn = 0;
			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
			{									// Finding NR in the room
				if (p->length1 > 1.5 * nrlength)continue;
				flag1 = findseq(nrseq, nrlength, p);
				if (flag1 == 0 && p->length2 == 0) nr_turn += NR_TURN;
			}
			np_bef = np_arr[1][y][x];
			np_arr[1][y][x] = 0;
			for (k = 0; k < np_bef; k++)
			{
				rntsyn = randd();
				flagntsyn = 1;
				if (rntsyn < PNF) flagntsyn = 0; // Uncatalyzed
				else if (nr_turn > 0)
				{
					nr_turn--;
					if (rntsyn < PNFR) flagntsyn = 0; // Catalyzed
				}
				if (flagntsyn == 0)
				{				      // Forming a nucleotide 
					p3 = (struct rna*)malloc(LEN);
					if (!p3) { printf("\t%d form_monomer--memeout\n", k + 1); exit(0); }
					randnt = randl(4) + 1;
					switch (randnt)
					{
					case 1:  p3->information[0][0] = A; break;
					case 2:  p3->information[0][0] = C; break;
					case 3:  p3->information[0][0] = G; break;
					case 4:  p3->information[0][0] = U; break;
					default: printf("form randnt error");
					}
					p3->information[0][1] = 0;
					p3->information[1][0] = 0;
					p3->length1 = 1;
					p3->length2 = 0;
					p3->prior = room_head[0][y][x];
					p3->next = room_head[0][y][x]->next;
					if (p3->next != room_head[0][y][x])(p3->next)->prior = p3;
					room_head[0][y][x]->next = p3;
				}
				else np_arr[0][y][x]++;
			}
		}
	}
	//======================= End of "Events of molecules"

	//======================================== Formation and break of protocells
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			if (m_arr[y][x] == 0)    // No membrane
			{					 // Possibly (cell) forming	
				if (fa_arr[0][y][x] + 2 * phl_arr[0][y][x] >= LAM)
				{
					if (randd() < 1 - pow(1 - PMF, fa_arr[0][y][x] + 2 * phl_arr[0][y][x] - LAM + 1))
					{
						m_arr[y][x] = fa_arr[0][y][x] + 2 * phl_arr[0][y][x];
						phl_on_m_arr[y][x] = phl_arr[0][y][x];
						fa_arr[0][y][x] = 0;
						phl_arr[0][y][x] = 0;
					}
				}
			}
			else     // Membrane existing
			{		 // Possibly (cell) breaking	
				if (m_arr[y][x] < LAM || randd() < PCB)
				{
					phl_arr[0][y][x] += phl_on_m_arr[y][x];
					fa_arr[0][y][x] += m_arr[y][x] - 2 * phl_on_m_arr[y][x];
					phl_on_m_arr[y][x] = 0;
					m_arr[y][x] = 0;
				}
			}
		}
	}
	//========================= End of "Formation and break of protocells"

	//==================================== Movement of protocells (including cell fusion)
	for (y = 0; y < N; y++)   // Initially flagging all rooms as "having not been considered"
		for (x = 0; x < N; x++)
			flagcmov[y][x] = 0;
	avail_xy_init();      // Initialization for xy_choose
	for (d = 0; d < ROOMNUM; d++)
	{
		xy_choose();    // Picks a room at random 
		if (m_arr[y][x] > 0 && flagcmov[y][x] == 0 && randd() < PMC)
		{
			randcase = randl(4);   // Four possible directions
			switch (randcase)
			{
			case 0:  // To left  
				if (m_arr[y][(N + x - 1) % N] > 0) // There is a protocell on the left
				{						  // Possibly fusing to Left
					if (randd() < PCF)
					{
						cell_fusing(0, -1);
						flagcmov[y][(N + x - 1) % N] = 1;
					}
				}
				else    // There is not a protocell on the left
				{		// Possibly moving to left
					left = 1; up = 1; down = 1;
					if (m_arr[y][(N + x - 2) % N] > 0) left = 0;
					if (m_arr[(N + y - 1) % N][(N + x - 1) % N] > 0) up = 0;
					if (m_arr[(y + 1) % N][(N + x - 1) % N] > 0) down = 0;
					if (left == 0 && up == 0 && down == 0)break; // No way for outer fatty acid moving, so no cell moving
					if (left == 1 && up == 0 && down == 0) oneway_outer_moving(0, 0, -1, -1);  //  Outer material moving allowed only to left
					if (left == 0 && up == 1 && down == 0) oneway_outer_moving(0, -1, -1, 0);  // Only up
					if (left == 0 && up == 0 && down == 1) oneway_outer_moving(0, 1, -1, 0);   // Only down
					if (left == 1 && up == 1 && down == 0) twoway_outer_moving(0, 0, -1, -1, -1, 0); // Left and up
					if (left == 1 && up == 0 && down == 1) twoway_outer_moving(0, 0, 1, -1, -1, 0);  // Left and down
					if (left == 0 && up == 1 && down == 1) twoway_outer_moving(0, -1, 1, -1, 0, 0);  // Up and down
					if (left == 1 && up == 1 && down == 1) threeway_outer_moving(0, 0, -1, 1, -1, -1, 0, 0); // Left, up and down
					cell_moving(0, -1); // Moving to left
					flagcmov[y][(N + x - 1) % N] = 1; // Flagging this room as "having been considered already"
				}
				break;
			case 1:    // To right
				if (m_arr[y][(x + 1) % N] > 0)
				{
					if (randd() < PCF)
					{
						cell_fusing(0, 1);
						flagcmov[y][(x + 1) % N] = 1;
					}
				}
				else
				{
					right = 1; up = 1; down = 1;
					if (m_arr[y][(x + 2) % N] > 0) right = 0;
					if (m_arr[(N + y - 1) % N][(x + 1) % N] > 0) up = 0;
					if (m_arr[(y + 1) % N][(x + 1) % N] > 0) down = 0;
					if (right == 0 && up == 0 && down == 0)break;
					if (right == 1 && up == 0 && down == 0) oneway_outer_moving(0, 0, 1, 1);
					if (right == 0 && up == 1 && down == 0) oneway_outer_moving(0, -1, 1, 0);
					if (right == 0 && up == 0 && down == 1) oneway_outer_moving(0, 1, 1, 0);
					if (right == 1 && up == 1 && down == 0) twoway_outer_moving(0, 0, -1, 1, 1, 0);
					if (right == 1 && up == 0 && down == 1) twoway_outer_moving(0, 0, 1, 1, 1, 0);
					if (right == 0 && up == 1 && down == 1) twoway_outer_moving(0, -1, 1, 1, 0, 0);
					if (right == 1 && up == 1 && down == 1) threeway_outer_moving(0, 0, -1, 1, 1, 1, 0, 0);
					cell_moving(0, 1);
					flagcmov[y][(x + 1) % N] = 1;
				}
				break;
			case 2:    // To up
				if (m_arr[(N + y - 1) % N][x] > 0)
				{
					if (randd() < PCF)
					{
						cell_fusing(-1, 0);
						flagcmov[(N + y - 1) % N][x] = 1;
					}
				}
				else
				{
					up = 1; left = 1; right = 1;
					if (m_arr[(N + y - 2) % N][x] > 0) up = 0;
					if (m_arr[(N + y - 1) % N][(N + x - 1) % N] > 0) left = 0;
					if (m_arr[(N + y - 1) % N][(x + 1) % N] > 0) right = 0;
					if (up == 0 && left == 0 && right == 0)break;
					if (up == 1 && left == 0 && right == 0) oneway_outer_moving(-1, -1, 0, 0);
					if (up == 0 && left == 1 && right == 0) oneway_outer_moving(-1, 0, 0, -1);
					if (up == 0 && left == 0 && right == 1) oneway_outer_moving(-1, 0, 0, 1);
					if (up == 1 && left == 1 && right == 0) twoway_outer_moving(-1, -1, 0, 0, 0, -1);
					if (up == 1 && left == 0 && right == 1) twoway_outer_moving(-1, -1, 0, 0, 0, 1);
					if (up == 0 && left == 1 && right == 1) twoway_outer_moving(-1, 0, 0, 0, -1, 1);
					if (up == 1 && left == 1 && right == 1) threeway_outer_moving(-1, -1, 0, 0, 0, 0, -1, 1);
					cell_moving(-1, 0);
					flagcmov[(N + y - 1) % N][x] = 1;
				}
				break;
			case 3:    // To down
				if (m_arr[(y + 1) % N][x] > 0)
				{
					if (randd() < PCF)
					{
						cell_fusing(1, 0);
						flagcmov[(y + 1) % N][x] = 1;
					}
				}
				else
				{
					down = 1; left = 1; right = 1;
					if (m_arr[(y + 2) % N][x] > 0) down = 0;
					if (m_arr[(y + 1) % N][(N + x - 1) % N] > 0) left = 0;
					if (m_arr[(y + 1) % N][(x + 1) % N] > 0) right = 0;
					if (down == 0 && left == 0 && right == 0)break;
					if (down == 1 && left == 0 && right == 0) oneway_outer_moving(1, 1, 0, 0);
					if (down == 0 && left == 1 && right == 0) oneway_outer_moving(1, 0, 0, -1);
					if (down == 0 && left == 0 && right == 1) oneway_outer_moving(1, 0, 0, 1);
					if (down == 1 && left == 1 && right == 0) twoway_outer_moving(1, 1, 0, 0, 0, -1);
					if (down == 1 && left == 0 && right == 1) twoway_outer_moving(1, 1, 0, 0, 0, 1);
					if (down == 0 && left == 1 && right == 1) twoway_outer_moving(1, 0, 0, 0, -1, 1);
					if (down == 1 && left == 1 && right == 1) threeway_outer_moving(1, 1, 0, 0, 0, 0, -1, 1);
					cell_moving(1, 0);
					flagcmov[(y + 1) % N][x] = 1;
				}
				break;
			default: printf("cell moving error");
			}
		}
	}
	//========================= End of "Movement of protocells (including cell fusion)"

	//========================================== Division of protocells
	for (y = 0; y < N; y++)   // Initially flagging all rooms as "having not been considered"
		for (x = 0; x < N; x++)
			flagcdiv[y][x] = 0;
	avail_xy_init();      // Initialization for xy_choose
	for (d = 0; d < ROOMNUM; d++)
	{
		xy_choose();    // Picks a room at random 
		if (m_arr[y][x] > 0 && flagcdiv[y][x] == 0 && randd() < PCD * (1 - 2.0 * LAM / m_arr[y][x]))
		{
			cell_div_times++;
			randcase = randl(4);   // Four possible directions
			switch (randcase)
			{
			case 0:    // To left
				if (m_arr[y][(N + x - 1) % N] > 0) break;
				left = 1; up = 1; down = 1;
				if (m_arr[y][(N + x - 2) % N] > 0) left = 0;
				if (m_arr[(N + y - 1) % N][(N + x - 1) % N] > 0) up = 0;
				if (m_arr[(y + 1) % N][(N + x - 1) % N] > 0) down = 0;
				if (left == 0 && up == 0 && down == 0)break; // No way for outer fatty acid moving, so no cell division
				if (left == 1 && up == 0 && down == 0) oneway_outer_moving(0, 0, -1, -1); // Outer material moving allowed only to left
				if (left == 0 && up == 1 && down == 0) oneway_outer_moving(0, -1, -1, 0); // Only up
				if (left == 0 && up == 0 && down == 1) oneway_outer_moving(0, 1, -1, 0);  // Only down
				if (left == 1 && up == 1 && down == 0) twoway_outer_moving(0, 0, -1, -1, -1, 0); // Left and up
				if (left == 1 && up == 0 && down == 1) twoway_outer_moving(0, 0, 1, -1, -1, 0);  // Left and down
				if (left == 0 && up == 1 && down == 1) twoway_outer_moving(0, -1, 1, -1, 0, 0);  // Up and down
				if (left == 1 && up == 1 && down == 1) threeway_outer_moving(0, 0, -1, 1, -1, -1, 0, 0); // Left, up and down
				cell_dividing(0, -1);      // Deviding to left
				flagcdiv[y][(N + x - 1) % N] = 1; // Flagging this room as "having been considered already"
				break;
			case 1: // To right
				if (m_arr[y][(x + 1) % N] > 0) break;
				right = 1; up = 1; down = 1;
				if (m_arr[y][(x + 2) % N] > 0) right = 0;
				if (m_arr[(N + y - 1) % N][(x + 1) % N] > 0) up = 0;
				if (m_arr[(y + 1) % N][(x + 1) % N] > 0) down = 0;
				if (right == 0 && up == 0 && down == 0)break;
				if (right == 1 && up == 0 && down == 0) oneway_outer_moving(0, 0, 1, 1);
				if (right == 0 && up == 1 && down == 0) oneway_outer_moving(0, -1, 1, 0);
				if (right == 0 && up == 0 && down == 1) oneway_outer_moving(0, 1, 1, 0);
				if (right == 1 && up == 1 && down == 0) twoway_outer_moving(0, 0, -1, 1, 1, 0);
				if (right == 1 && up == 0 && down == 1) twoway_outer_moving(0, 0, 1, 1, 1, 0);
				if (right == 0 && up == 1 && down == 1) twoway_outer_moving(0, -1, 1, 1, 0, 0);
				if (right == 1 && up == 1 && down == 1) threeway_outer_moving(0, 0, -1, 1, 1, 1, 0, 0);
				cell_dividing(0, 1);
				flagcdiv[y][(x + 1) % N] = 1;
				break;
			case 2: // To up
				if (m_arr[(N + y - 1) % N][x] > 0) break;
				up = 1; left = 1; right = 1;
				if (m_arr[(N + y - 2) % N][x] > 0) up = 0;
				if (m_arr[(N + y - 1) % N][(N + x - 1) % N] > 0) left = 0;
				if (m_arr[(N + y - 1) % N][(x + 1) % N] > 0) right = 0;
				if (up == 0 && left == 0 && right == 0)break;
				if (up == 1 && left == 0 && right == 0) oneway_outer_moving(-1, -1, 0, 0);
				if (up == 0 && left == 1 && right == 0) oneway_outer_moving(-1, 0, 0, -1);
				if (up == 0 && left == 0 && right == 1) oneway_outer_moving(-1, 0, 0, 1);
				if (up == 1 && left == 1 && right == 0) twoway_outer_moving(-1, -1, 0, 0, 0, -1);
				if (up == 1 && left == 0 && right == 1) twoway_outer_moving(-1, -1, 0, 0, 0, 1);
				if (up == 0 && left == 1 && right == 1) twoway_outer_moving(-1, 0, 0, 0, -1, 1);
				if (up == 1 && left == 1 && right == 1) threeway_outer_moving(-1, -1, 0, 0, 0, 0, -1, 1);
				cell_dividing(-1, 0);
				flagcdiv[(N + y - 1) % N][x] = 1;
				break;
			case 3:    // To down
				if (m_arr[(y + 1) % N][x] > 0) break;
				down = 1; left = 1; right = 1;
				if (m_arr[(y + 2) % N][x] > 0) down = 0;
				if (m_arr[(y + 1) % N][(N + x - 1) % N] > 0) left = 0;
				if (m_arr[(y + 1) % N][(x + 1) % N] > 0) right = 0;
				if (down == 0 && left == 0 && right == 0)break;
				if (down == 1 && left == 0 && right == 0) oneway_outer_moving(1, 1, 0, 0);
				if (down == 0 && left == 1 && right == 0) oneway_outer_moving(1, 0, 0, -1);
				if (down == 0 && left == 0 && right == 1) oneway_outer_moving(1, 0, 0, 1);
				if (down == 1 && left == 1 && right == 0) twoway_outer_moving(1, 1, 0, 0, 0, -1);
				if (down == 1 && left == 0 && right == 1) twoway_outer_moving(1, 1, 0, 0, 0, 1);
				if (down == 0 && left == 1 && right == 1) twoway_outer_moving(1, 0, 0, 0, -1, 1);
				if (down == 1 && left == 1 && right == 1) threeway_outer_moving(1, 1, 0, 0, 0, 0, -1, 1);
				cell_dividing(1, 0);
				flagcdiv[(y + 1) % N][x] = 1;
				break;
			default: printf("cell division error");
			}
		}
	}
	//=========================== End of "Division of protocells"
}

/////////////////////////////////////////////////////////////////////////
void record(void)        // Data recording at every interval step (RECINT) 
{
	int ch_num[MAX_RNA_LENGTH], ch_nr_num[MAX_RNA_LENGTH], ch_gr_num[MAX_RNA_LENGTH], ch_nrgr_num[MAX_RNA_LENGTH], ch_0_num[MAX_RNA_LENGTH], over_twenty, si;

	FILE* fptxt, * fptxt1;
	errno_t err, err1;
	err = fopen_s(&fptxt, "file.txt", "at");
	if (err != 0) { printf("cannot open file");  exit(-1); }
	err1 = fopen_s(&fptxt1, "file1.txt", "at");
	if (err1 != 0) { printf("cannot open file1");  exit(-1); }

	int	t_phl_on_m, t_m, t_phl_on_m_g, t_m_g, t_phl_on_m_ng, t_m_ng, gr_count;
	int t_np_in, t_np_in_nr, t_np_in_no_nr, t_npp_in, t_npp_in_nr, t_npp_in_no_nr, t_gp_in, t_np_out, t_npp_out, t_gp_out;
	float c_g_num, c_ng_num, c_num, c_nr_num, c_no_nr_num, no_c_num;

	gr_num[g] = 0;         // Number of GR
	nr_num[g] = 0;         // Number of NR
	npr_num[g] = 0;        // Number of NPR
	contr1_num[g] = 0;     // Number of RNAs containing control-1
	contr2_num[g] = 0;
	gr_incell_num[g] = 0;  // Number of GR in protocells
	nr_incell_num[g] = 0;  // Number of NR in protocells
	npr_incell_num[g] = 0; // Number of NPR in protocells
	contr1_incell_num[g] = 0;  // Number of RNA containing control-1 in protocells
	contr2_incell_num[g] = 0;

	total_nt_mat[g] = 0;   // Total materials for nucleotide-precursor's precursors, nucleotide precursors, nucleotides and RNAs (quotients in measurement of nucleotides)
	total_fa_mat[g] = 0;   // Total materials for fatty acids and fatty acid residues in phospholipids (quotients in measurement of fatty acids)
	total_g_mat[g] = 0;    // Total materials for glycerophosphate precursors, glycerophosphates and the glycerophosphate portion in phospholipids (quotients in measurement of glycerophosphate)

	npp_num[g] = 0;        // Number of nucleotide-precursor's precursors
	np_num[g] = 0;         // Number of nucleotide precursors
	rna_num[g] = 0;        // Number of nucleotides and RNAs
	fa_num[g] = 0;         // Number of fatty acids out of the membrane

	gp_num[g] = 0;         // Number of glycerophosphate precursors 
	g_num[g] = 0;          // Number of glycerophosphates
	phl_num[g] = 0;        // Number of phospholipids out of the membrane
	phl_on_m_num[g] = 0;   // Number of phospholipids within the membrane

	cell_num[g] = 0;       // Number of protocells
	cell_0_num[g] = 0;     // Number of protocells without functional RNAs
	cell_gr_num[g] = 0;	   // Number of protocells containing GR
	cell_nr_num[g] = 0;	   // Number of protocells containing NR
	cell_npr_num[g] = 0;   // Number of protocells containing NPR
	cell_nrgr_num[g] = 0;  // Number of protocells containing both NR and GR
	cell_nrnpr_num[g] = 0; // Number of protocells containing both NR and NPR
	cell_grnpr_num[g] = 0; // Number of protocells containing both GR and NPR
	cell_nrgrnpr_num[g] = 0; // Number of protocells containing all NR, GR and NPR

	cell_ctrl_num[g] = 0;  // Number of protocells containing control-1

	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			np_num[g] += np_arr[0][y][x];
			npp_num[g] += npp_arr[0][y][x];
			fa_num[g] += fa_arr[0][y][x];

			gp_num[g] += gp_arr[0][y][x];
			g_num[g] += g_arr[0][y][x];
			phl_num[g] += phl_arr[0][y][x];
			phl_on_m_num[g] += phl_on_m_arr[y][x];

			total_fa_mat[g] += fa_arr[0][y][x] + m_arr[y][x] + 2 * phl_arr[0][y][x];

			total_g_mat[g] += gp_arr[0][y][x] + g_arr[0][y][x] + phl_arr[0][y][x] + phl_on_m_arr[y][x];

			flagnr = 1; flagnpr = 1; flaggr = 1;  flagctrl = 1;
			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
			{
				rna_num[g]++;
				total_nt_mat[g] += p->length1 + p->length2;
				flag1 = findseq(nrseq, nrlength, p);
				if (flag1 == 0)
				{
					nr_num[g]++;
					flagnr = 0;
					if (m_arr[y][x] > 0) nr_incell_num[g]++;
				}
				flag2 = findseq(nprseq, nprlength, p);
				if (flag2 == 0)
				{
					npr_num[g]++;
					flagnpr = 0;
					if (m_arr[y][x] > 0) npr_incell_num[g]++;
				}
				flag3 = findseq(contrseq1, contrlength1, p);
				if (flag3 == 0)
				{
					contr1_num[g]++;
					flagctrl = 0;
					if (m_arr[y][x] > 0) contr1_incell_num[g]++;
				}
				flag3 = findseq(contrseq2, contrlength2, p);
				if (flag3 == 0)
				{
					contr2_num[g]++;
					if (m_arr[y][x] > 0) contr2_incell_num[g]++;
				}
				flag3 = findseq(grseq, grlength, p);
				if (flag3 == 0)
				{
					gr_num[g]++;
					flaggr = 0;
					if (m_arr[y][x] > 0) gr_incell_num[g]++;
				}
			}
			if (m_arr[y][x] > 0)
			{
				cell_num[g]++;
				if (flagnr != 0 && flaggr != 0 && flagnpr != 0)cell_0_num[g]++;
				if (flagnr == 0 && flaggr != 0 && flagnpr != 0)cell_nr_num[g]++;
				if (flagnr != 0 && flaggr == 0 && flagnpr != 0)cell_gr_num[g]++;
				if (flagnr != 0 && flaggr != 0 && flagnpr == 0)cell_npr_num[g]++;
				if (flagnr == 0 && flaggr == 0 && flagnpr != 0)cell_nrgr_num[g]++;
				if (flagnr == 0 && flaggr != 0 && flagnpr == 0)cell_nrnpr_num[g]++;
				if (flagnr != 0 && flaggr == 0 && flagnpr == 0)cell_grnpr_num[g]++;
				if (flagnr == 0 && flaggr == 0 && flagnpr == 0)cell_nrgrnpr_num[g]++;
				if (flagctrl == 0)cell_ctrl_num[g]++;
			}
		}
	}
	total_nt_mat[g] += np_num[g] + npp_num[g];

	printf("\n%d: nr=(%d)%d, gr=(%d)%d, contr1=(%d)%d, contr2=(%d)%d, npr=(%d)%d\n[tn=%d,rna=%d,np=%d,npp=%d,ta=%d,a=%d,tg=%d,g=%d,gp=%d,phl=%d(%d),c=%d,c0=%d,cnr=%d,cgr=%d,cnpr=%d,cnrgr=%d,cnrnpr=%d,cgrnpr=%d,cnrgrnpr=%d,cctrl=%d]\n",
		i, (int)nr_incell_num[g], (int)nr_num[g], (int)gr_incell_num[g], (int)gr_num[g], (int)contr1_incell_num[g], (int)contr1_num[g], (int)contr2_incell_num[g], (int)contr2_num[g],
		(int)npr_incell_num[g], (int)npr_num[g],
		(int)total_nt_mat[g], (int)rna_num[g], (int)np_num[g], (int)npp_num[g], (int)total_fa_mat[g], (int)fa_num[g],
		(int)total_g_mat[g], (int)g_num[g], (int)gp_num[g], (int)phl_num[g], (int)phl_on_m_num[g],
		(int)cell_num[g], (int)cell_0_num[g], (int)cell_nr_num[g], (int)cell_gr_num[g], (int)cell_npr_num[g], (int)cell_nrgr_num[g], (int)cell_nrnpr_num[g], (int)cell_grnpr_num[g], (int)cell_nrgrnpr_num[g], (int)cell_ctrl_num[g]);

	if (i == 0)fprintf(fptxt, "\n\n");
	fprintf(fptxt, "\n%d: nr=(%d)%d, gr=(%d)%d, contr1=(%d)%d, contr2=(%d)%d, npr=(%d)%d\n[tn=%d,rna=%d,np=%d,npp=%d,ta=%d,a=%d,tg=%d,g=%d,gp=%d,phl=%d(%d),c=%d,c0=%d,cnr=%d,cgr=%d,cnpr=%d,cnrgr=%d,cnrnpr=%d,cgrnpr=%d,cnrgrnpr=%d,cctrl=%d]\n",
		i, (int)nr_incell_num[g], (int)nr_num[g], (int)gr_incell_num[g], (int)gr_num[g], (int)contr1_incell_num[g], (int)contr1_num[g], (int)contr2_incell_num[g], (int)contr2_num[g],
		(int)npr_incell_num[g], (int)npr_num[g],
		(int)total_nt_mat[g], (int)rna_num[g], (int)np_num[g], (int)npp_num[g], (int)total_fa_mat[g], (int)fa_num[g],
		(int)total_g_mat[g], (int)g_num[g], (int)gp_num[g], (int)phl_num[g], (int)phl_on_m_num[g],
		(int)cell_num[g], (int)cell_0_num[g], (int)cell_nr_num[g], (int)cell_gr_num[g], (int)cell_npr_num[g], (int)cell_nrgr_num[g], (int)cell_nrnpr_num[g], (int)cell_grnpr_num[g], (int)cell_nrgrnpr_num[g], (int)cell_ctrl_num[g]);

	fprintf(fptxt1, "step=%d: nr=%d,  gr=%d, contr1=%d, contr2=%d, npr=%d, c=%d, c0=%d,cnr=%d,cgr=%d,cnpr=%d,cnrgr=%d,cnrnpr=%d,cgrnpr=%d,cnrgrnpr=%d,cctrl=%d, \n",
		i, (int)nr_num[g], (int)gr_num[g], (int)contr1_num[g], (int)contr2_num[g], (int)npr_num[g],
		(int)cell_num[g], (int)cell_0_num[g], (int)cell_nr_num[g], (int)cell_gr_num[g], (int)cell_npr_num[g], (int)cell_nrgr_num[g], (int)cell_nrnpr_num[g], (int)cell_grnpr_num[g], (int)cell_nrgrnpr_num[g], (int)cell_ctrl_num[g]);

	for (si = 0; si < MAX_RNA_LENGTH; si++)
	{
		ch_num[si] = 0;
		ch_nr_num[si] = 0;
		ch_gr_num[si] = 0;
		ch_nrgr_num[si] = 0;
		ch_0_num[si] = 0;
	}
	over_twenty = 0;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
			{
				ch_num[p->length1 - 1]++;
				if (p->length1 > 20)
				{
					over_twenty++;
					for (int t = 0; t < p->length1; t++)   // Reading out the RNA sequence longer than 20 nt
					{
						switch (p->information[0][t])
						{
						case 1: printf("A"); fprintf(fptxt, "A"); break;
						case 2: printf("C"); fprintf(fptxt, "C"); break;
						case 3: printf("G"); fprintf(fptxt, "G"); break;
						case 4: printf("U"); fprintf(fptxt, "U"); break;
						default: printf("error");
						}
					}
					printf("\t%d\n", p->length1); fprintf(fptxt, "\t%d\n", p->length1);
					for (int t = 0; t < p->length2; t++)
					{
						switch (p->information[1][t])
						{
						case 1: printf("A"); fprintf(fptxt, "A"); break;
						case 2: printf("C"); fprintf(fptxt, "C"); break;
						case 3: printf("G"); fprintf(fptxt, "G"); break;
						case 4: printf("U"); fprintf(fptxt, "U"); break;
						default: printf("error");
						}
					}
					printf("\t%d\n\n", p->length2); fprintf(fptxt, "\t%d\n\n", p->length2);
				}
				flag1 = findseq(nrseq, nrlength, p);
				flag3 = findseq(grseq, grlength, p);
				if (flag1 == 0 && flag3 != 0)ch_nr_num[p->length1 - 1]++;
				else if (flag1 != 0 && flag3 == 0)ch_gr_num[p->length1 - 1]++;
				else if (flag1 == 0 && flag3 == 0)ch_nrgr_num[p->length1 - 1]++;
				else ch_0_num[p->length1 - 1]++;
			}
		}
	}

	for (si = 0; si < 20; si++)   // Show relevant information for RNA sequences no longer than 20 nt
	{
		printf("%d(%d|%d/%d^%d), ", ch_num[si], ch_nr_num[si], ch_gr_num[si], ch_nrgr_num[si], ch_0_num[si]);
		fprintf(fptxt, "%d(%d|%d/%d^%d), ", ch_num[si], ch_nr_num[si], ch_gr_num[si], ch_nrgr_num[si], ch_0_num[si]);
	}
	printf("\nover20=%d   step=%d\n", over_twenty, i);	fprintf(fptxt, "\nover20=%d   step=%d\n", over_twenty, i);

	t_phl_on_m = 0;   // For summarizing detailed information about protocells and molecules
	t_m = 0;

	c_ng_num = 0;
	t_phl_on_m_ng = 0;
	t_m_ng = 0;
	c_g_num = 0;
	t_phl_on_m_g = 0;
	t_m_g = 0;

	t_np_in = 0;
	t_np_in_nr = 0;
	t_np_in_no_nr = 0;
	t_npp_in = 0;
	t_npp_in_nr = 0;
	t_npp_in_no_nr = 0;
	t_gp_in = 0;
	t_np_out = 0;
	t_npp_out = 0;
	t_gp_out = 0;
	c_num = 0;
	c_nr_num = 0;
	c_no_nr_num = 0;
	no_c_num = 0;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			if (m_arr[y][x] == 0)   // Summarizing detailed information about protocells and molecules
			{
				no_c_num++;
				t_np_out += np_arr[0][y][x];
				t_npp_out += npp_arr[0][y][x];
				t_gp_out += gp_arr[0][y][x];
			}
			else
			{
				c_num++;
				t_np_in += np_arr[0][y][x];
				t_npp_in += npp_arr[0][y][x];
				t_gp_in += gp_arr[0][y][x];

				flagnr = 1;
				for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
				{
					flag1 = findseq(nrseq, nrlength, p);
					if (flag1 == 0) { flagnr = 0; break; }
				}
				if (flagnr == 0)
				{
					c_nr_num++;
					t_np_in_nr += np_arr[0][y][x];
					t_npp_in_nr += npp_arr[0][y][x];
				}
				else
				{
					c_no_nr_num++;
					t_np_in_no_nr += np_arr[0][y][x];
					t_npp_in_no_nr += npp_arr[0][y][x];
				}
			}
			if (m_arr[y][x] == 0) fprintf(fptxt, "0\t");  // Showing relevant information concerning the spatial distribution
			else
			{
				gr_count = 0;
				for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
				{
					if (findseq(grseq, grlength, p) == 0)gr_count++;
				}
				fprintf(fptxt, "%d/%d-%d\t", phl_on_m_arr[y][x], m_arr[y][x], gr_count);

				t_phl_on_m += phl_on_m_arr[y][x];
				t_m += m_arr[y][x];
				if (gr_count != 0)
				{
					c_g_num++;
					t_phl_on_m_g += phl_on_m_arr[y][x];
					t_m_g += m_arr[y][x];
				}
				else
				{
					c_ng_num++;
					t_phl_on_m_ng += phl_on_m_arr[y][x];
					t_m_ng += m_arr[y][x];
				}
			}
		}
		fprintf(fptxt, "\n");
	}
	if (t_m == 0) { printf("no cell\n"); fprintf(fptxt, "no cell\n"); }   // Show the information about protocells and molecules summarized above
	else
	{
		printf("2*phl/m=%f, 2*phl/m-no_gr=%f, 2*phl/m-gr=%f, c_size=%f, c_size_ng=%f, c_size_g=%f, cell_div_times=%ld\n",
			2.0 * t_phl_on_m / t_m, 2.0 * t_phl_on_m_ng / t_m_ng, 2.0 * t_phl_on_m_g / t_m_g,
			t_m / (c_g_num + c_ng_num), t_m_ng / c_ng_num, t_m_g / c_g_num, cell_div_times);
		fprintf(fptxt, "2*phl/m=%f, 2*phl/m-no_gr=%f, 2*phl/m-gr=%f, c_size=%f, c_size_ng=%f, c_size_g=%f, cell_div_times=%ld\n",
			2.0 * t_phl_on_m / t_m, 2.0 * t_phl_on_m_ng / t_m_ng, 2.0 * t_phl_on_m_g / t_m_g,
			t_m / (c_g_num + c_ng_num), t_m_ng / c_ng_num, t_m_g / c_g_num, cell_div_times);
	}
	printf("-----c_num=%f, c_nr_num=%f, c_no_nr_num=%f, t_np_in=%d, t_np_in_nr=%d, t_np_in_no_nr=%d, t_npp_in=%d, t_gp_in=%d, no_c_num=%f, t_np_out=%d, t_npp_out=%d, t_gp_out=%d\n",
		c_num, c_nr_num, c_no_nr_num, t_np_in, t_np_in_nr, t_np_in_no_nr, t_npp_in, t_gp_in, no_c_num, t_np_out, t_npp_out, t_gp_out);
	fprintf(fptxt, "-----c_num=%f, c_nr_num=%f, c_no_nr_num=%f, t_np_in=%d, t_np_in_nr=%d, t_np_in_no_nr=%d, t_npp_in=%d, t_gp_in=%d, no_c_num=%f, t_np_out=%d, t_npp_out=%d, t_gp_out=%d\n",
		c_num, c_nr_num, c_no_nr_num, t_np_in, t_np_in_nr, t_np_in_no_nr, t_npp_in, t_gp_in, no_c_num, t_np_out, t_npp_out, t_gp_out);
	if (c_num == 0)
	{
		printf("-----Concentration: np=(N/A)%f, npp=(N/A)%f, gp=(N/A)%f\n", t_np_out / no_c_num, t_npp_out / no_c_num, t_gp_out / no_c_num);
		fprintf(fptxt, "-----Concentration: np=(N/A)%f, npp=(N/A)%f, gp=(N/A)%f\n", t_np_out / no_c_num, t_npp_out / no_c_num, t_gp_out / no_c_num);
	}
	else
	{
		printf("-----Concentration: np=(%f-%f/%f)%f, npp=(%f-%f/%f)%f, gp=(%f)%f\n",
			t_np_in / c_num, t_np_in_nr / c_nr_num, t_np_in_no_nr / c_no_nr_num, t_np_out / no_c_num, t_npp_in / c_num, t_npp_in_nr / c_nr_num, t_npp_in_no_nr / c_no_nr_num, t_npp_out / no_c_num, t_gp_in / c_num, t_gp_out / no_c_num);
		fprintf(fptxt, "-----Concentration: np=(%f-%f/%f)%f, npp=(%f-%f/%f)%f, gp=(%f)%f\n",
			t_np_in / c_num, t_np_in_nr / c_nr_num, t_np_in_no_nr / c_no_nr_num, t_np_out / no_c_num, t_npp_in / c_num, t_npp_in_nr / c_nr_num, t_npp_in_no_nr / c_no_nr_num, t_npp_out / no_c_num, t_gp_in / c_num, t_gp_out / no_c_num);
	}

	g++;
	fclose(fptxt);
	fclose(fptxt1);
}


/////////////////////////////////////////////////////////////////////////
void freepool(void)        // Memory releasing  
{
	int m;
	for (m = 0; m < 2; m++)
	{
		for (y = 0; y < N; y++)
		{
			for (x = 0; x < N; x++)
			{
				while (1)
				{
					if (room_head[m][y][x]->next != room_head[m][y][x])
					{
						p = room_head[m][y][x]->next;
						room_head[m][y][x]->next = p->next;
						free(p);
					}
					else break;
				}
				free(room_head[m][y][x]);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////
int main()
{
	inits();            // Initialization of the system
	for (i = 0; i <= STEPNUM; i++)  // Monte-Carlo cycle
	{
		if (i == INOCUSTEP0) inocu_cell();  // Inoculation of empty protocells
		if (i == INOCUSTEP1) inoculate1();  // Inoculation of GR and Control-1 into empty protocells

		if (i >= STAREC && i % RECINT == 0)    // Data recording at every interval step (RECINT) 
		{
			record();
		}
		unit_action();		// Action of units (molecules and protocells) in the system
	}
	freepool();      // Memory releasing

	return (0);
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> End of the program







