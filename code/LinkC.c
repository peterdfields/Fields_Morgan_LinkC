/* Calculates linkage desequilibria between a couple of loci in several
 * subdivided populations, and Ohtas partitioning of variance of the
 * desequilibrium within and between sub-populations.
 *
 * Program transcribed manually from Pascal to C by
 * Mark Morgan (mmm2a@virginia.edu) and Peter Fields (pdf8z@virginia.edu). Questions concerning
 * LinkC should be directed to Peter Fields.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "LinkCConfig.h"
#include "LinkCTypes.h"
#include "LinkCUtilities.h"

#define STRING_SIZE 20
#define MAX_FILENAME_SIZE 2056

static matpop popn, nmax;
static matal alef;
static ptald alefld, homold, sumalf, homosm;
static ptdl sumld, totsum;
static matind ind;
static int npopct, indld, fixchk, nodat, sumind, negfre;
static matnumal nall, totn, ndr1, ndr2;
static ptdes burrow, corrad;
static float bursq, bursd, chisq, z, prob;
static float fone, fthree, fseven, seven;
static matdl chisc, rcprob;
static ptst wta;
static int npop, nloc;
static FILE *file1, *file2, *file3, *filet;
static char response;
static int pop, loc, loco, cc;
static int allchk, df;
static char **aloc;
static char **apop;
static char tempo[MAX_FILENAME_SIZE];

static float gauss(float x);
static float chiprb(float chisq, int df);
static void init(LinkCConfiguration config);
static void initcouple(LinkCConfiguration config, int loc1, int loc2);
static void lecture(LinkCConfiguration config);
static void remp(int lc1, int lc2);
static void ldcalc(int allchk, int l1, int l2);
static void clear(int a, int b);
static void trans(int c, int d);
static void alefp(LinkCConfiguration config);
static void freq(const char *outfreq);
static void checkprob(LinkCConfiguration config, int loci, int locii,
		int *choix);
static void ldlwrt(int loci, int locii, int allchk);
static void outpop(LinkCConfiguration config, int allchk, float signi, int loci,
		int locii);
static void ohta(int loc1, int loc2);
static void condldwrt(int lc1, int lc2);
static void ohtawrt();
static void linC(LinkCConfiguration config);

#define CONST_ZERO 0.000001f

float gauss(float x) {
	double z, y, w;

	if (!(x == CONST_ZERO)) {
		y = fabs(x) / 2;
		if (y >= 3)
			z = 1.0f;
		else if (y > 1) {
			y -= 2;
			z = -0.000045155659 * y + 0.000152529290;
			z = z * y - 0.000019538132;
			z = z * y - 0.000676904986;
			z = z * y + 0.001390604284;
			z = z * y - 0.000794620820;
			z = z * y - 0.002034254874;
			z = z * y + 0.006549791214;
			z = z * y - 0.010557625006;
			z = z * y + 0.011630447319;
			z = z * y - 0.009279453341;
			z = (z * y + 0.005353579108) * y - 0.002141268741;
			z = (z * y + 0.000535310849) * y + 0.999936657524;

			/* Was commented out in the original Pascal
			 z = (((((((((((((-0.045155659 * y
			 +0.152529290) * y - 0.019538132) * y
			 -0.676904986) * y + 1.390604284) * y
			 -0.794620820) * y - 2.034254874) * y
			 +6.549791214) * y - 10.557625006) * y
			 +11.630447319) * y - 9.279453341) * y
			 +5.353579108) * y - 2.141268741) * y
			 +0.535310849) * y + 999.936657524;
			 z /= 1000;
			 */
		} else {
			w = y * y;
			z = 0.000124818987 * w - 0.001075204047;
			z = z * w + 0.005198775019;
			z = (z * w - 0.019198292004) * w + 0.059054035642;
			z = (z * w - 0.151968751364) * w + 0.319152932694;
			z = ((z * w - 0.531923007300) * w + 0.797884560593) * y * 2;

			/* Was commented out in the original Pascal
			 z = ((((((((0.000124818987 * w
			 -0.001075204047) * w + 0.005198775019) * w
			 -0.019198292004) * w + 0.059054035642) * w
			 -0.151968751364) * w + 0.319152932694) * w
			 -0.531923007300) * w + 0.797884560593) * y * 2;
			 */

		}
	} else
		z = 0;

	if (x > CONST_ZERO)
		return (1 + z) / 2;
	else
		return (1 - z) / 2;
}

float chiprb(float chisq, int df) {
	int even, bigx;
	float a, s, z, e, c;
	float chisqb, y;
	float ret;

	if (chisq <= CONST_ZERO)
		return 0;

	if (df > 30)
		return 1
				- gauss(
						-1 * sqrt(4.5 * df)
								* (pow(chisq / df, 0.333333333) + 2 / (9 * df)
										- 1));

	a = 0.5 * chisq;
	ret = 0;
	even = ((df % 2) == 0);
	bigx = (chisq >= 200);
	if (df <= 2) {
		if (bigx)
			ret = 1;
		else {
			if (even)
				ret = 1 - exp(-a);
			else
				ret = 1 - 2 * gauss(-1 * sqrt(chisq));
		}
	} else {
		chisqb = 0.5 * (df - 1);
		if (bigx) {
			if (even) {
				s = 0;
				z = 0;
				e = 0;
			} else {
				s = 2 * gauss(-1 * sqrt(chisq));
				z = -0.5;
				e = 0.57236494295; /* log(sqrt(pi)) */
			}

			if (z <= 0)
				z = 0.000001;
			c = log(z) + e;
			z += 1;
			while (z < chisqb) {
				e = log(z + e);
				s = exp(c * z - a - e) + s;
				z = z + 1;
			}
			ret = 1 - s;
		} else {
			y = exp(-1 * a);
			if (even) {
				z = 0;
				s = y;
				e = 1;
			} else {
				s = 2 * gauss(-1 * sqrt(chisq));
				z = -0.5;
				e = 0.564189583548 / sqrt(a); /* 1 / sqrt(pi) */
			}

			c = 0;
			z = z + 1;
			while (z < chisqb) {
				e = e * a / z;
				c += e;
				z += 1;
			}
			ret = 1 - (c * y + s);
		}
	}

	return ret;
}

void init(LinkCConfiguration config) {
	int all;
	int popi, loci;

	memset(popn, 0, sizeof(int) * config.npopm);
	memset(nmax, 0, sizeof(int) * config.npopm);

	for (popi = 0; popi < config.npopm; popi++) {
		memset(ind[popi], 0, sizeof(int) * config.nlocm);

		for (loci = 0; loci < config.nlocm; loci++) {
			for (all = 0; all < config.alm; all++) {
				alef[popi][loci][all] = 0;
				(*wta)[loci][all] = 0;
			}
		}
	}

	memset(nall, 0, sizeof(int) * config.nlocm);
	memset(totn, 0, sizeof(int) * config.nlocm);
}

void initcouple(LinkCConfiguration config, int loc1, int loc2) {
	int al1;
	int al2;

	npopct = npop;
	indld = 0;
	fixchk = 0;
	nodat = 0;
	sumind = 0;
	negfre = 0;
	bursq = 0;
	bursd = 0;
	chisq = 0;
	z = 0;
	prob = 0;
	fone = 0;
	fthree = 0;
	fseven = 0;
	seven = 0;

	for (al1 = 0; al1 < config.alm; al1++) {
		for (al2 = 0; al2 < config.alm; al2++) {
			(*sumld)[al1][al2] = 0;
			(*totsum)[al1][al2] = 0;
			(*burrow)[al1][al2] = 0;
			(*corrad)[al1][al2] = 0;
			(*alefld)[al1][loc2] = 0;
			(*alefld)[al2][loc1] = 0;
			(*homold)[al1][loc2] = 0;
			(*homold)[al2][loc1] = 0;
			(*sumalf)[al1][loc2] = 0;
			(*sumalf)[al2][loc1] = 0;
			(*homosm)[al1][loc2] = 0;
			(*homosm)[al2][loc1] = 0;
			chisc[al1][al2] = 0;
			rcprob[al1][al2] = 0;
		}
	}
}

void lecture(LinkCConfiguration config) {
	int indi;
	int i;
	int j;
	int loc;

	for (indi = 1; indi <= nmax[pop - 1]; indi++) {
		for (loc = 1; loc <= nloc; loc++) {
			fscanf(file1, "%d %d", &ndr1[loc - 1], &ndr2[loc - 1]);

			for (i = 1; i <= nall[loc - 1]; i++) {
				if (ndr1[loc - 1] != 0) {
					if (ndr1[loc - 1] == i) {
						for (j = 1; j <= nall[loc - 1]; j++) {
							if (ndr2[loc - 1] != 0) {
								if (ndr2[loc - 1] == j) {
									if (i != j) {
										alef[pop - 1][loc - 1][i - 1] = 1
												+ alef[pop - 1][loc - 1][i - 1];

										alef[pop - 1][loc - 1][j - 1] = 1
												+ alef[pop - 1][loc - 1][j - 1];

										ind[pop - 1][loc - 1] = 1
												+ ind[pop - 1][loc - 1];
									} else if (i == j) {
										alef[pop - 1][loc - 1][i - 1] = 2
												+ alef[pop - 1][loc - 1][i - 1];

										ind[pop - 1][loc - 1] = 1
												+ ind[pop - 1][loc - 1];
									} /*cas o÷ i<>j ou i=j*/
								} /*condition ndr2[loc]=j*/
							} /*ndr2[loc]<>0*/
						} /*boucle sur ndr2[loc]*/
					} /*condition ndr1[loc]=i*/
				} /*condition ndr1[loc]<>0*/
			} /*boucle sur ndr1[loc]*/
		} /*boucle sur les locus*/
	} /*boucle sur les individus intra sous-pop*/
}

/**********************************************************************/
/****Remplissage des matrices                                **********/
/**********************************************************************/

void remp(int lc1, int lc2) {
	int i, j, k, l;
	int indi;

	/*Comptage des­rents types de gamÉtes par couple de locus, etc...*/
	for (indi = 1; indi <= nmax[pop - 1]; indi++) {
		for (i = 1; i <= nloc; i++)
			fscanf(file1, "%d %d", &(ndr1[i - 1]), &(ndr2[i - 1]));

		fscanf(file1, "\n");

		for (i = 1; i <= nall[lc1 - 1]; i++) {
			/*compteur 1er allÉle 1er locus*/
			if (ndr1[lc1 - 1] != 0) {
				if (ndr1[lc1 - 1] == i) {
					for (j = 1; j <= nall[lc1 - 1]; j++) {
						/*compteur 2Éme allÉle  1er locus*/
						if (ndr2[lc1 - 1] != 0) {
							if (ndr2[lc1 - 1] == j) {
								for (k = 1; k <= nall[lc2 - 1]; k++) {
									/*compteur 1er allÉle 2Éme locus*/
									if (ndr1[lc2 - 1] != 0) {
										if (ndr1[lc2 - 1] == k) {
											for (l = 1; l <= nall[lc2 - 1];
													l++) {
												/* compteur 2Éme allÉle 2Éme
												 * loc
												 */
												if (ndr2[lc2 - 1] != 0) {
													if (ndr2[lc2 - 1] == l) {
														/* Unindented for readability */
														/*cas o÷ i=j et k=l*/
														if ((i == j)
																&& (k == l)) {
															(*sumld)[i - 1][k
																	- 1] =
																	2
																			+ ((*sumld)[i
																					- 1][k
																					- 1]);
															(*alefld)[i - 1][lc2
																	- 1] =
																	2
																			+ ((*alefld)[i
																					- 1][lc2
																					- 1]);
															(*alefld)[k - 1][lc1
																	- 1] =
																	2
																			+ ((*alefld)[k
																					- 1][lc1
																					- 1]);
															(*homold)[i - 1][lc2
																	- 1] =
																	1
																			+ ((*homold)[i
																					- 1][lc2
																					- 1]);
															(*homold)[k - 1][lc1
																	- 1] =
																	1
																			+ ((*homold)[k
																					- 1][lc1
																					- 1]);
															indld = 1 + indld;
															(*totsum)[i - 1][k
																	- 1] =
																	2
																			+ ((*totsum)[i
																					- 1][k
																					- 1]);
															(*sumalf)[i - 1][lc2
																	- 1] =
																	2
																			+ ((*sumalf)[i
																					- 1][lc2
																					- 1]);
															(*sumalf)[k - 1][lc1
																	- 1] =
																	2
																			+ ((*sumalf)[k
																					- 1][lc1
																					- 1]);
															(*homosm)[i - 1][lc2
																	- 1] =
																	1
																			+ ((*homosm)[i
																					- 1][lc2
																					- 1]);
															(*homosm)[k - 1][lc1
																	- 1] =
																	1
																			+ ((*homosm)[k
																					- 1][lc1
																					- 1]);
															sumind = 1 + sumind;
														}

														/*cas o÷ i=j et k<>l*/
														if ((i == j)
																&& (k != l)) {
															(*sumld)[i - 1][k
																	- 1] =
																	1
																			+ ((*sumld)[i
																					- 1][k
																					- 1]);
															(*sumld)[i - 1][l
																	- 1] =
																	1
																			+ ((*sumld)[i
																					- 1][l
																					- 1]);
															(*alefld)[i - 1][lc2
																	- 1] =
																	2
																			+ ((*alefld)[i
																					- 1][lc2
																					- 1]);
															(*alefld)[k - 1][lc1
																	- 1] =
																	1
																			+ ((*alefld)[k
																					- 1][lc1
																					- 1]);
															(*alefld)[l - 1][lc1
																	- 1] =
																	1
																			+ ((*alefld)[l
																					- 1][lc1
																					- 1]);
															(*homold)[i - 1][lc2
																	- 1] =
																	1
																			+ ((*homold)[i
																					- 1][lc2
																					- 1]);
															indld = 1 + indld;
															(*totsum)[i - 1][k
																	- 1] =
																	1
																			+ ((*totsum)[i
																					- 1][k
																					- 1]);
															(*totsum)[i - 1][l
																	- 1] =
																	1
																			+ ((*totsum)[i
																					- 1][l
																					- 1]);

															(*sumalf)[i - 1][lc2
																	- 1] =
																	2
																			+ ((*sumalf)[i
																					- 1][lc2
																					- 1]);
															(*sumalf)[k - 1][lc1
																	- 1] =
																	1
																			+ ((*sumalf)[k
																					- 1][lc1
																					- 1]);
															(*sumalf)[l - 1][lc1
																	- 1] =
																	1
																			+ ((*sumalf)[l
																					- 1][lc1
																					- 1]);
															(*homosm)[i - 1][lc2
																	- 1] =
																	1
																			+ ((*homosm)[i
																					- 1][lc2
																					- 1]);
															sumind = 1 + sumind;
														}

														/*cas o÷ i<>j et k=l*/
														if ((i != j)
																&& (k == l)) {
															(*sumld)[i - 1][k
																	- 1] =
																	1
																			+ ((*sumld)[i
																					- 1][k
																					- 1]);
															(*sumld)[j - 1][k
																	- 1] =
																	1
																			+ ((*sumld)[j
																					- 1][k
																					- 1]);

															(*alefld)[i - 1][lc2
																	- 1] =
																	1
																			+ ((*alefld)[i
																					- 1][lc2
																					- 1]);
															(*alefld)[j - 1][lc2
																	- 1] =
																	1
																			+ ((*alefld)[j
																					- 1][lc2
																					- 1]);
															(*alefld)[k - 1][lc1
																	- 1] =
																	2
																			+ ((*alefld)[k
																					- 1][lc1
																					- 1]);
															(*homold)[k - 1][lc1
																	- 1] =
																	1
																			+ ((*homold)[k
																					- 1][lc1
																					- 1]);
															indld = 1 + indld;
															(*totsum)[i - 1][k
																	- 1] =
																	1
																			+ ((*totsum)[i
																					- 1][k
																					- 1]);
															(*totsum)[j - 1][l
																	- 1] =
																	1
																			+ ((*totsum)[j
																					- 1][l
																					- 1]);

															(*sumalf)[i - 1][lc2
																	- 1] =
																	1
																			+ ((*sumalf)[i
																					- 1][lc2
																					- 1]);
															(*sumalf)[j - 1][lc2
																	- 1] =
																	1
																			+ ((*sumalf)[j
																					- 1][lc2
																					- 1]);
															(*sumalf)[k - 1][lc1
																	- 1] =
																	2
																			+ ((*sumalf)[k
																					- 1][lc1
																					- 1]);
															(*homosm)[k - 1][lc1
																	- 1] =
																	1
																			+ ((*homosm)[k
																					- 1][lc1
																					- 1]);
															sumind = 1 + sumind;
														}

														/*cas o÷ i<>j et k<>l*/
														if ((i != j)
																&& (k != l)) {
															/*on doit d­clarer sumld en*/
															(*sumld)[i - 1][k
																	- 1] =
																	0.5
																			+ (*sumld)[i
																					- 1][k
																					- 1];

															/*integer*/
															(*sumld)[j - 1][k
																	- 1] =
																	0.5
																			+ (*sumld)[j
																					- 1][k
																					- 1];
															(*sumld)[i - 1][l
																	- 1] =
																	0.5
																			+ (*sumld)[i
																					- 1][l
																					- 1];
															(*sumld)[j - 1][l
																	- 1] =
																	0.5
																			+ (*sumld)[j
																					- 1][l
																					- 1];
															(*alefld)[i - 1][lc2
																	- 1] =
																	1
																			+ (*alefld)[i
																					- 1][lc2
																					- 1];
															(*alefld)[j - 1][lc2
																	- 1] =
																	1
																			+ (*alefld)[j
																					- 1][lc2
																					- 1];
															(*alefld)[k - 1][lc1
																	- 1] =
																	1
																			+ (*alefld)[k
																					- 1][lc1
																					- 1];
															(*alefld)[l - 1][lc1
																	- 1] =
																	1
																			+ (*alefld)[l
																					- 1][lc1
																					- 1];
															indld = 1 + indld;

															/*mÃme rem pour totsum que*/
															(*totsum)[i - 1][k
																	- 1] =
																	0.5
																			+ (*totsum)[i
																					- 1][k
																					- 1];

															/*pour sumld*/
															(*totsum)[j - 1][k
																	- 1] =
																	0.5
																			+ (*totsum)[j
																					- 1][k
																					- 1];
															(*totsum)[i - 1][l
																	- 1] =
																	0.5
																			+ (*totsum)[i
																					- 1][l
																					- 1];

															(*totsum)[j - 1][l
																	- 1] =
																	0.5
																			+ (*totsum)[j
																					- 1][l
																					- 1];

															(*sumalf)[i - 1][lc2
																	- 1] =
																	1
																			+ (*sumalf)[i
																					- 1][lc2
																					- 1];
															(*sumalf)[j - 1][lc2
																	- 1] =
																	1
																			+ (*sumalf)[j
																					- 1][lc2
																					- 1];
															(*sumalf)[k - 1][lc1
																	- 1] =
																	1
																			+ (*sumalf)[k
																					- 1][lc1
																					- 1];
															(*sumalf)[l - 1][lc1
																	- 1] =
																	1
																			+ (*sumalf)[l
																					- 1][lc1
																					- 1];
															sumind = 1 + sumind;
														}

														/* Indent them original ones back or something ? */

													} /*condition ndr2[lc2]=l*/
												} /*condition ndr2[lc2]<>0*/
											} /* coptur 2Éme allÉle 2Éme locus*/
										} /* condition ndr1[lc2]=k*/
									} /*condition ndr1[lc2]<>0*/
								} /*compteur 1er allÉle 2Éme locus*/
							} /*condition ndr2[lc1]=j*/
						} /*condition ndr2[lc1]<>0*/
					} /*compteur 2Éme allÉle 1er locus*/
				} /*condition ndr1[lc1]=i*/
			} /*condition ndr1[lc1]<>0*/
		} /*compteur 1er allÉle 1er locus*/
	} /*boucle sur les individus*/

	/* boucle sur les sous-pop et les couples de
	 * locus  dans le program principal
	 */
} /* remp */

/**********************************************************************/
/**********   Calcul du desequilibre de liaison pop/pop   *************/
/**********************************************************************/
void ldcalc(int allchk, int l1, int l2) {
	int a1, a2;
	int ngfrc, ddl1, ddl2, d1, d2;
	float sumbsq, sumde, sumz;
	float zstat, bsq, chi, zbar;
	float alefa, assoc, p, q, alefb;
	int bool;
	float hwcora, hwcorb, denom;

	sumbsq = 0;
	sumde = 0;
	sumz = 0;
	ngfrc = 0;
	denom = 0;

	if (indld != 0) {
		bool = 0;

		for (a1 = 1; a1 <= nall[l1 - 1]; a1++) {
			alefa = (float) ((*alefld)[a1 - 1][l2 - 1]) / (2 * indld);

			if ((alefa != 0)) {
				if ((alefa == 1)) {
					bool = 1;

					for (a2 = 1; a2 <= nall[l2 - 1]; a2++) {
						assoc = (float) ((*alefld)[a2 - 1][l1 - 1]) / indld;
						fone = fone + powf((float) assoc, (float) 2) * indld;
						alefa = (float) ((*alefld)[a2 - 1][l1 - 1])
								/ (2 * indld);
						fseven = fseven + powf(alefa, (float) 2) * indld;

						if (allchk == 0) {
							p = alef[pop - 1][l2 - 1][a2 - 1];

							seven = (p * p) + seven;
							fthree = fthree + assoc * indld * alefa;
						}
					}
				} else /*alefa=1*/
				{ /*alefa < 1*/
					for (a2 = 1; a2 <= nall[l2 - 1]; a2++) {
						alefb = (float) ((*alefld)[a2 - 1][l1 - 1])
								/ (2 * indld);

						if (alefb != 0) { /*aucun des allÉles fixe*/
							if (alefb < 1) {
								assoc = (float) ((*sumld)[a1 - 1][a2 - 1])
										/ indld;
								fone = fone + powf(assoc, (float) 2) * indld;
								fthree = fthree + assoc * indld * alefa * alefb;
								fseven = fseven
										+ powf(alefa, (float) 2)
												* powf(alefb, (float) 2)
												* indld;

								if (allchk == 0) {
									p = alef[pop - 1][l1 - 1][a1 - 1];
									q = alef[pop - 1][l2 - 1][a2 - 1];
									seven = (p * p) * (q * q) + seven;
								}

								hwcora = ((float) (*homold)[a1 - 1][l2 - 1]
										/ indld) - powf(alefa, (float) 2);
								hwcorb = ((float) (*homold)[a2 - 1][l1 - 1]
										/ indld) - powf(alefb, (float) 2);

								(*burrow)[a1 - 1][a2 - 1] = (assoc
										- 2 * alefa * alefb)
										* ((float) indld / (indld - 1));

								denom = (alefa * (1 - alefa) + hwcora)
										* (alefb * (1 - alefb) + hwcorb);

								if (denom != 0) {
									(*corrad)[a1 - 1][a2 - 1] =
											(float) ((*burrow)[a1 - 1][a2 - 1])
													/ sqrt(denom);

									if ((*corrad)[a1 - 1][a2 - 1] >= 1)
										(*corrad)[a1 - 1][a2 - 1] = 0.99999;

									if ((*corrad)[a1 - 1][a2 - 1] < -1)
										(*corrad)[a1 - 1][a2 - 1] = -0.99999;

									zstat =
											(float) (0.5
													* log(
															(1
																	+ (*corrad)[a1
																			- 1][a2
																			- 1])
																	/ (1
																			- (*corrad)[a1
																					- 1][a2
																					- 1])));
									zstat = zstat * (indld - 3);

									sumz = sumz + fabs(zstat);

									bsq = fabs((*burrow)[a1 - 1][a2 - 1])
											* fabs((*burrow)[a1 - 1][a2 - 1]);

									chi = (float) (indld * bsq)
											/ (alefa * alefb);
									ngfrc += 1;
									sumbsq += bsq;
									chisq += chi;
									sumde += denom;
								} /* denom<>0*/
							} else /*alefb <1*/
							{ /*alefb=1*/
								bool = 1;

								for (a1 = 1; a1 <= nall[l1 - 1]; a1++) {
									assoc = (float) ((*alefld)[a1 - 1][l2 - 1])
											/ indld;
									fone += powf(assoc, (float) 2) * indld;

									alefa = (float) ((*alefld)[a1 - 1][l2 - 1])
											/ (2 * indld);
									fseven += powf(alefa, (float) 2) * indld;

									if (allchk == 0) {
										p = alef[pop - 1][l1 - 1][a1 - 1];
										seven = (p * p) + seven;
										fthree = fthree + assoc * indld * alefa;
									}
								}
							} /*aleb=1*/
						} /*alefb<>0*/
					} /*boucle sur a2*/
				} /*alefa<1*/
			} /*alefa<>0*/
		} /*boucle sur a1*/

		if (bool == 0) {
			if (sumde != 0) {
				bursq = sqrt((float) sumbsq / ngfrc);
				bursd = (float) sumbsq / sumde;
				if ((indld - 3) > 0) {
					zbar = (float) sumz / ((indld - 3) * ngfrc);

					if (exp(2 * zbar) != 1) {
						z = (float) (exp(2 * zbar) - 1) / (exp(2 * zbar) + 1);

						/*calcul des ddl pour les tests de chi2*/
						ddl1 = 0;
						ddl2 = 0;

						for (d1 = 1; d1 <= nall[l1 - 1]; d1++) {
							if ((*alefld)[d1 - 1][l2 - 1] != 0)
								ddl1 += 1;
						}

						for (d2 = 1; d2 <= nall[l2 - 1]; d2++) {
							if ((*alefld)[(d2) - 1][l1 - 1] != 0)
								ddl2 += 1;
						}

						negfre = (ddl1 - 1) * (ddl2 - 1);
						df = negfre;

						/*def. du parametre df pour chiprb,outpop*/
						prob = fabs(1 - chiprb(chisq, df));

						if (prob <= 0.0001)
							prob = 0.0001;
					} /*condition exp(2*zbar<>1*/
				} /*condition indld<>3*/
			} /*condition  sumde<>0*/
		} else /*bool=false*/
		{
			fixchk = 1;
		} /*condition 1 allÉle fix­ · au moins un des 2 locus*/

	} else /*condition indld[l1,l2]<>0*/
	{ /*un des locus est manquant*/
		nodat = 1;
		npopct = npopct - 1;
	}
} /*ldcalc*/

/******************************************************************/
/******R­initialisations des matrices pour ldcalc pop/pop *********/
/******************************************************************/
void clear(int a, int b) {
	int aa, bb;

	indld = 0;
	fixchk = 0;
	nodat = 0;
	negfre = 0;
	bursq = 0;
	bursd = 0;
	chisq = 0;
	z = 0;
	prob = 0;

	for (aa = 1; aa <= nall[a - 1]; aa++) {
		for (bb = 1; bb <= nall[b - 1]; bb++) {
			(*sumld)[aa - 1][bb - 1] = 0;
			(*burrow)[aa - 1][bb - 1] = 0;
			(*corrad)[aa - 1][bb - 1] = 0;
			(*alefld)[aa - 1][b - 1] = 0;
			(*alefld)[bb - 1][a - 1] = 0;
			(*homold)[aa - 1][b - 1] = 0;
			(*homold)[bb - 1][a - 1] = 0;

			chisc[aa - 1][bb - 1] = 0;
			rcprob[aa - 1][bb - 1] = 0;
		}
	}
} /*clear*/

/******************************************************************/
/*transformation des matrices pour ldcalc dans toute la population*/
/******************************************************************/
void trans(int c, int d) {
	int cc, dd;

	indld = sumind;

	for (cc = 1; cc <= nall[c - 1]; cc++) {
		for (dd = 1; dd <= nall[d - 1]; dd++) {
			(*sumld)[cc - 1][dd - 1] = (*totsum)[cc - 1][dd - 1];
			(*alefld)[cc - 1][d - 1] = (*sumalf)[cc - 1][d - 1];
			(*alefld)[dd - 1][c - 1] = (*sumalf)[dd - 1][c - 1];
			(*homold)[cc - 1][d - 1] = (*homosm)[cc - 1][d - 1];
			(*homold)[dd - 1][c - 1] = (*homosm)[dd - 1][c - 1];
		}
	}
} /*trans*/

/******************************************************************/
/**********Calcul des fr­quences all­liques pop/pop****************/
/******************************************************************/
void alefp(LinkCConfiguration config) {
	ptstat wgt;
	int pp, ll, ic, ip;
	int ipp, icc, ial;
	int jpp, jcc, jal;

	/* This is no longer sufficient, we also have to create the matrix
	 * itself.
	 */
	wgt = (ptstat) malloc(1 * sizeof(matstat)); /*init¿ de la matrice wgt*/
	(*wgt) = createMatstat(config);

	for (pp = 1; pp <= npop; pp++) {
		for (ll = 1; ll <= nloc; ll++) {
			(*wgt)[pp - 1][(ll) - 1] = 0.0;
		} /*calcul du nb total d'individus dans les populations*/
	}

	for (ic = 1; ic <= nloc; ic++) {
		for (ip = 1; ip <= npop; ip++)
			totn[ic - 1] += ind[ip - 1][ic - 1];
	} /*calcul des freq. all­liques et idem pond­r­es*/

	for (ipp = 1; ipp <= npop; ipp++) {
		for (icc = 1; icc <= nloc; icc++) {
			for (ial = 1; ial <= nall[icc - 1]; ial++) {
				if (ind[ipp - 1][icc - 1] != 0)
					alef[ipp - 1][icc - 1][ial - 1] /= 2
							* (ind[ipp - 1][icc - 1]);
			}
		}
	}

	for (jcc = 1; jcc <= nloc; jcc++) {
		for (jal = 1; jal <= nall[jcc - 1]; jal++) {
			for (jpp = 1; jpp <= npop; jpp++) {
				if (ind[jpp - 1][jcc - 1] != 0) {
					(*wgt)[jpp - 1][jcc - 1] = (float) (ind[jpp - 1][jcc - 1])
							/ (totn[jcc - 1]);

					(*wta)[jcc - 1][jal - 1] = ((*wta)[jcc - 1][jal - 1])
							+ (alef[jpp - 1][jcc - 1][jal - 1]
									* ((*wgt)[jpp - 1][jcc - 1]));
				}
			}
		}
	}

	destroyMatstat(config, *wgt);
	free(wgt);
}

/******************************************************************/
/***********Edition du tableau des fr­quences all­liques***********/
/******************************************************************/
void freq(const char *outfreq) {
	FILE *file4;
	int al, m;

	file4 = fopen(outfreq, "w");

	if (file4 == NULL ) {
		fprintf(stderr, "Unable to open file \"%s\" for writing.\n", outfreq);
		exit(-1);
	}

	fprintf(file4, "\n");
	fprintf(file4, "ALLELIC FREQUENCIES IN SUBPOPULATIONS\n");
	fprintf(file4, "\n");
	fprintf(file4, "\n");
	fprintf(file4, "LOCUS  TOT");

	for (m = 1; m <= 6; m++) {
		pop = (10 * m - 9);

		if (pop <= npop) {
			while ((pop <= (10 * m)) && (pop <= npop)) {
				fprintf(file4, "    %2d  ", pop);
				pop = pop + 1;
			}

			fprintf(file4, "\n");
			for (loc = 1; loc <= nloc; loc++) {
				fprintf(file4, "%s\n", aloc[loc - 1]);
				fprintf(file4, " (N)  %4d", totn[loc - 1]);
				pop = (10 * m - 9);
				while ((pop <= (10 * m)) && (pop <= npop)) {
					fprintf(file4, "    %3d", ind[pop - 1][loc - 1]);

					pop += 1;
				}

				fprintf(file4, "\n");
				for (al = 1; al <= nall[loc - 1]; al++) {
					fprintf(file4, "%3d  %5.3f", al, ((*wta)[loc - 1][al - 1]));
					pop = (10 * m - 9);

					while ((pop <= (10 * m)) && (pop <= npop)) {
						fprintf(file4, "  %5.3f",
								alef[pop - 1][loc - 1][al - 1]);

						pop += 1;
					} /*boucle sur les populations*/

					fprintf(file4, "\n");

				} /*boucle sur les alleles*/

				fprintf(file4, "\n");

			} /*boucle sur les locus*/

			fprintf(file4, "\n");
		} /*condition sur pop*/
	} /*boucle sur m*/

	fclose(file4);
} /*freq*/

/******************************************************************/
/*********    Decision si                                **********/
/***Ecriture detaillee des LD par couple d'alleles quand        ***/
/**** -le test global pour un couple de locus est significatif*****/
/**** -un couple d"alleles est fortement correle              *****/
/******************************************************************/
void checkprob(LinkCConfiguration config, int loci, int locii, int *choix) {
	int al1, al2, count;

	count = 0;
	if (nodat == 0) {
		for (al1 = 1; al1 <= nall[loci - 1]; al1++) {
			for (al2 = 1; al2 <= nall[locii - 1]; al2++) {
				if ((*burrow)[al1 - 1][al2 - 1] != 0) {
					chisc[al1 - 1][al2 - 1] = (fabs(
							powf(((*corrad)[al1 - 1][al2 - 1]), 2))) * indld;

					df = 1;
					rcprob[al1 - 1][al2 - 1] = (fabs(
							1 - chiprb(chisc[al1 - 1][al2 - 1], df)));

					if ((rcprob[al1 - 1][al2 - 1] <= 0.0001))
						rcprob[al1 - 1][al2 - 1] = 0.0001;

					if ((rcprob[al1 - 1][al2 - 1] <= config.signi))
						count += 1;
				}
			}
		}

		if (count >= 1)
			*choix = 1;
		else
			*choix = 0;
	}
} /*checkprob*/

/******************************************************************/
/***Ecriture detaillee des LD par couple d'alleles quand        ***/
/**** -le test global pour un couple de locus est significatif*****/
/**** -un couple d"alleles est fortement correle              *****/
/******************************************************************/
void ldlwrt(int loci, int locii, int allchk) {
	int bb, al1, al2;

	fprintf(filet, "\n");
	fprintf(filet, "\n");
	for (bb = 1; bb <= 80; bb++)
		fprintf(filet, "*");

	fprintf(filet, "\n");
	fprintf(filet, "detailed output %s - %s\n", aloc[loci - 1],
			aloc[locii - 1]);
	fprintf(filet, "\n");

	if ((allchk == 0) && (nodat == 0)) {
		fprintf(filet, "population %s\n", apop[pop - 1]);

		for (bb = 1; bb <= 80; bb++)
			fprintf(filet, "*");

		fprintf(filet, "\n");
		fprintf(filet, "\n");
		fprintf(filet, "Sample size= %2d individuals\n", indld);
		fprintf(filet, "\n");
		fprintf(filet, "   ALLELES         T(IJ)      D(IJ)        R(IJ)  ");
		fprintf(filet, " CHI-SQUARE   PROB.\n");

		for (al1 = 1; al1 <= nall[loci - 1]; al1++) {
			for (al2 = 1; al2 <= nall[locii - 1]; al2++) {
				if (((*burrow)[al1 - 1][al2 - 1]) != 0) {
					fprintf(filet, "  %2d  -  %2d        ", al1, al2);
					fprintf(filet, "%5.1f    ", (*sumld)[al1 - 1][al2 - 1]);
					fprintf(filet, "%8.5f    ", (*burrow)[al1 - 1][al2 - 1]);
					fprintf(filet, "%8.5f    %7.2f",
							(*corrad)[al1 - 1][al2 - 1],
							chisc[al1 - 1][al2 - 1]);

					fprintf(filet, "      %6.4f\n", rcprob[al1 - 1][al2 - 1]);
				}
			} /*boucle sur al2*/
		} /*boucle sur al1*/

		fprintf(filet, "\n");
		fprintf(filet, "variance     = %7.5f\n", bursd);
		fprintf(filet, "common correlation   = %7.5f\n", z);
		fprintf(filet, "chi-square(%d)       = %4.2f\n", negfre, chisq);
		fprintf(filet, "prob.            = %6.4f\n", prob);
	} /*condition allchk=0 et nodat=0*/

	if ((allchk == 1) && (nodat == 0)) {
		fprintf(filet, "all subpopulations\n");

		for (bb = 1; bb <= 80; bb++)
			fprintf(filet, "*");

		fprintf(filet, "\n");
		fprintf(filet, "\n");
		fprintf(filet, "Sample size= %2d individuals\n", indld);
		fprintf(filet, "\n");
		fprintf(filet, "   ALLELES         T(IJ)      D(IJ)        R(IJ)  ");
		fprintf(filet, " CHI-SQUARE   PROB.\n");

		for (al1 = 1; al1 <= nall[loci - 1]; al1++) {
			for (al2 = 1; al2 <= nall[locii - 1]; al2++) {
				if (((*burrow)[al1 - 1][al2 - 1]) != 0) {
					fprintf(filet, "  %2d  -  %2d        ", al1, al2);
					fprintf(filet, "%5.1f    ", (*sumld)[al1 - 1][al2 - 1]);
					fprintf(filet, "%8.5f    ", (*burrow)[al1 - 1][al2 - 1]);
					fprintf(filet, "%8.5f    %7.2f",
							(*corrad)[al1 - 1][al2 - 1],
							chisc[al1 - 1][al2 - 1]);
					fprintf(filet, "      %6.4f\n", rcprob[al1 - 1][al2 - 1]);
				}
			} /*boucle sur al2*/
		} /*boucle sur al1*/

		fprintf(filet, "\n");
		fprintf(filet, "variance     = %7.5f\n", bursd);
		fprintf(filet, "common correlation   = %7.5f\n", z);
		fprintf(filet, "chi-square(%2d)      = %4.2f\n", negfre, chisq);
		fprintf(filet, "prob.            = %6.4f\n", prob);
	} /*condition allchk=1 et nodat=0*/
}

/******************************************************************/
/***************Ecriture condensee des LD**************************/
/******************************************************************/

void outpop(LinkCConfiguration config, int allchk, float signi, int loci,
		int locii) {
	int sig;
	int choose;

	if (fixchk >= 1) {
		if (allchk == 0)
			fprintf(file2, "%-4s             ", apop[pop - 1]);
		else
			fprintf(file2, "ALL SUBPOP       ");

		fprintf(file2, ": AN ALLELE AT ONE OR BOTH LOCI IS FIXED\n");
	} else if (nodat == 1) {
		if (allchk == 0)
			fprintf(file2, "%-4s             ", apop[pop - 1]);
		else
			fprintf(file2, "ALL SUBPOP       ");

		fprintf(file2, "NO INFORMATION WAS COLLECTED UPON ONE OF THE LOCI\n");
	} else {
		if (allchk == 0)
			fprintf(file2, "%-4s             ", apop[pop - 1]);
		else
			fprintf(file2, "ALL SUBPOP       ");

		fprintf(file2, "    %3d        ", indld);
		fprintf(file2, "  %7.5f      ", z);
		fprintf(file2, " %5.2f     ", chisq);
		fprintf(file2, "%3d     %6.4f\n", negfre, prob);
	}

	sig = 0;

	if (fixchk == 0) {
		if (nodat == 0) {
			checkprob(config, loci, locii, &choose);
			if ((prob <= signi) || (choose == 1)) {
				sig = 1;
				ldlwrt(loci, locii, allchk);
			}
		}
	}
}

/******************************************************************/
/*********Ecriture des composantes de la variance du LD************/
/******************************************************************/

void ohta(int loc1, int loc2) {
	float twelve, ftwo, fsix, fthirt, thirt;
	float dit, dstp, disp, dis, dst;
	float xbar, ybar, gbar;
	float fthre, fon, fsevn, sevn, diff;
	int al1, al2;

	if (npopct > 1) {
		twelve = 0;
		ftwo = 0;
		fsix = 0;
		fthirt = 0;
		thirt = 0;

		for (al1 = 1; al1 <= nall[loc1 - 1]; al1++) {
			for (al2 = 1; al2 <= nall[loc2 - 1]; al2++) {
				xbar = (float) (*sumalf)[al1 - 1][loc2 - 1] / (2 * sumind);
				ybar = (float) (*sumalf)[al2 - 1][loc1 - 1] / (2 * sumind);
				gbar = (float) (*totsum)[al1 - 1][al2 - 1] / sumind;

				ftwo += (gbar * gbar);
				fsix += 2 * xbar * ybar * gbar;
				fthirt += 4 * (xbar * xbar) * (ybar * ybar);
				thirt += 4
						* (((*wta)[loc1 - 1][al1 - 1])
								* ((*wta)[loc1 - 1][al1 - 1]))
						* (((*wta)[loc2 - 1][al2 - 1])
								* ((*wta)[loc2 - 1][al2 - 1]));

				for (pop = 1; pop <= npop; pop++) {
					if (((alef[pop - 1][loc1 - 1][al1 - 1] <= 1)
							&& (alef[pop - 1][loc2 - 1][al2 - 1] <= 1))) {
						diff = alef[pop - 1][loc1 - 1][al1 - 1]
								* alef[pop - 1][loc2 - 1][al2 - 1]
								* (*wta)[loc1 - 1][al1 - 1]
								* (*wta)[loc2 - 1][al2 - 1];

						twelve += diff;
					} /*alleles inf a 1*/
				} /*boucle populations*/
			} /*boucle alleles*/
		}

		fthre = (float) (2 * fthree) / sumind;
		fon = (float) fone / sumind;
		fsevn = (float) (4 * fseven) / sumind;
		sevn = (float) (4 * seven) / npopct;
		twelve = (float) (4 * twelve) / npopct;
		dit = fon + fthirt - 2 * fsix;
		dstp = ftwo + fthirt - 2 * fsix;
		disp = fon - ftwo;
		dis = fon + fsevn - 2 * fthre;
		dst = sevn + thirt - 2 * twelve;
		dst = (float) dst / 4;

		if ((dit <= 0) && (dit >= -0.000001))
			dit = 0;
		if ((dis <= 0) && (dis >= -0.000001))
			dis = 0;
		if ((disp <= 0) && (disp >= -0.000001))
			disp = 0;
		if ((dst <= 0) && (dst >= -0.000001))
			dst = 0;
		if ((dstp <= 0) && (dstp >= -0.000001))
			dstp = 0;

		fprintf(file3, "  %s - %s    ", aloc[loc1 - 1], aloc[loc2 - 1]);
		fprintf(file3, "    %7.5f  %7.5f     ", dis, disp);
		fprintf(file3, "   %7.5f  %7.5f      ", dst, dstp);
		fprintf(file3, "   %7.5f\n", dit);
	} /*npopct superieur a 1*/
}

/*****************************************************************/
/**********  En tetes des fichiers de sortie   *******************/
/*****************************************************************/

/**********Ecriture condensee des LD ***********/
void condldwrt(int lc1, int lc2) {
	int a;

	fprintf(file2, "\n");

	for (a = 1; a <= 80; a++)
		fprintf(file2, "*");

	fprintf(file2, "\n");
	fprintf(file2, "ANALYSIS OF LINKAGE DISEQUILIBRIUM BETWEEN %s - %s \n",
			aloc[lc1 - 1], aloc[lc2 - 1]);
	fprintf(file2, "\n");

	for (a = 1; a <= 80; a++)
		fprintf(file2, "*");

	fprintf(file2, "\n");
	fprintf(file2, "\n");
	fprintf(file2, "                 NUMBER OF      COMMON         CHI-\n");
	fprintf(file2, "  POPULATION      COMPARISONS    CORRELATION    SQUARE");
	fprintf(file2, "    D.F.    PROB.\n");

	for (a = 1; a <= 80; a++)
		fprintf(file2, "-");

	fprintf(file2, "\n");
}

/************** Ecriture des composantes de la variance du deseq ***/
void ohtawrt() {
	int aa;

	for (aa = 1; aa <= 80; aa++)
		fprintf(file3, "*");

	fprintf(file3, "\n");
	fprintf(file3, "VARIANCE COMPONENTS OF LINKAGE DISEQUILIBRIUM\n");

	for (aa = 1; aa <= 80; aa++)
		fprintf(file3, "-");

	fprintf(file3, "\n");
	fprintf(file3, "                   WITHIN SUBPOPULATION    ");
	fprintf(file3, "BETWEEN SUBPOPULATION    TOTAL POPULATION");
	fprintf(file3, "\n");
	fprintf(file3, "LOCI COMPARED         COMPONENTS              COMPONENT\n");
	fprintf(file3, "                   --------------------    ");
	fprintf(file3, "---------------------     ----------------\n");
	fprintf(file3,
			"                     D(IS)      D\"(IS)       D(ST)      D\"(ST)");
	fprintf(file3, "          D(IT)\n");

	for (aa = 1; aa <= 80; aa++)
		fprintf(file3, "-");

	fprintf(file3, "\n");
}

/******************************************************************/
/***********************Programme principal************************/
/******************************************************************/

void linC(LinkCConfiguration config) {
	char tmp[STRING_SIZE];

	strcpy(tempo, config.input);
	strcat(tempo, ".tmp");

	/* ------------ Malloc memory for pointers ------------- */
	alefld = (ptald) malloc(sizeof(matald));
	*alefld = createMatald(config);

	homold = (ptald) malloc(sizeof(matald));
	*homold = createMatald(config);

	sumalf = (ptald) malloc(sizeof(matald));
	*sumalf = createMatald(config);

	homosm = (ptald) malloc(sizeof(matald));
	*homosm = createMatald(config);

	sumld = (ptdl) malloc(sizeof(matdl));
	*sumld = createMatdl(config);

	totsum = (ptdl) malloc(sizeof(matdl));
	*totsum = createMatdl(config);

	burrow = (ptdes) malloc(sizeof(matdes));
	*burrow = createMatdes(config);

	corrad = (ptdes) malloc(sizeof(matdes));
	*corrad = createMatdes(config);

	wta = (ptst) malloc(sizeof(matst));
	*wta = createMatst(config);

	/* ----------- initialise pointers --------------------- */
	init(config);

	/* ----------- opening input file to read ----------------- */
	if ((file1 = fopen(config.input, "r")) != NULL ) {
		fscanf(file1, "%d%d\n", &npop, &nloc);

		/*lecture des noms de populations et du nb d'ind/population*/
		for (pop = 1; pop <= npop; pop++) {
			fgets(apop[pop - 1], 256, file1);
			apop[pop - 1][strlen(apop[pop - 1]) - 1] = '\0';

			fscanf(file1, "%d\n", &(nmax[pop - 1]));
		} /*lecture des noms des locus et du nb d'allÉles/locus*/

		for (loc = 1; loc <= nloc; loc++) {
			fscanf(file1, "%s %d\n", aloc[loc - 1], &(nall[loc - 1]));
		} /*lecture des g­notypes par pop,locus et couple de locus*/

		for (pop = 1; pop <= npop; pop++)
			lecture(config);

		alefp(config); /*calcul des frequences alleliques*/

		if (config.freqopt == 1)
			freq(config.outfreq);

		fclose(file1);
	} /* file 1 */

	if ((config.ldopt == 1) || (config.ldtotopt == 1)) {
		file2 = fopen(config.outdes, "w");
		if (!file2) {
			fprintf(stderr, "Unable to open output file \"%s\" for writing.\n",
					config.outdes);
			exit(-1);
		}
	} /*ouverture des fichiers*/

	if (config.ohtaopt == 1) {
		file3 = fopen(config.outohta, "w");
		if (!file3) {
			fprintf(stderr, "Unable to open output file \"%s\" for writing.\n",
					config.outohta);
			exit(-1);
		}

		ohtawrt();
	}

	for (loc = 1; loc <= (nloc - 1); loc++) {
		for (loco = (loc + 1); loco <= nloc; loco++) {
			/*boucle sur les couples de locus*/

			initcouple(config, loc, loco);
			if ((config.ldopt == 1) || (config.ldtotopt == 1)) {
				filet = fopen(tempo, "w");
				if (!filet) {
					fprintf(stderr, "Unable to open \"%s\" for writing.\n",
							config.outohta);
					exit(-1);
				}

				condldwrt(loc, loco);
			} /* end if */

			file1 = fopen(config.input, "r");
			if (!file1) {
				fprintf(stderr, "Unable to open \"%s\" for reading.\n",
						config.input);
				exit(-1);
			}

			for (cc = 1; cc <= (1 + (2 * npop) + nloc); cc++) {
				/* this is a dummy one .. no need */
				fgets(tmp, 256, file1);
			}

			for (pop = 1; pop <= npop; pop++) {
				allchk = 0;
				remp(loc, loco);

				if ((config.ldopt == 1) || (config.ohtaopt == 1)) {
					ldcalc(allchk, loc, loco);

					if (config.ldopt == 1)
						outpop(config, allchk, config.signi, loc, loco);
				} /* end if */

				clear(loc, loco); /*reinitialisation des matrices */
				/*pour le test du chi2*/
			} /* end for *//*boucle sur les pop*/

			fclose(file1);

			if (config.ohtaopt == 1) {
				allchk = 0;
				ohta(loc, loco);
			} /*ohta*/

			if (config.ldtotopt == 1) {
				allchk = 1;
				trans(loc, loco); /* transformation des matrices
				 indld en sumind et alefld en sumalf*/
				ldcalc(allchk, loc, loco);
				outpop(config, allchk, config.signi, loc, loco);
			} /* end if */

			if ((config.ldopt == 1) || (config.ldtotopt == 1)) {
				fclose(filet);
				filet = fopen(tempo, "r");
				if (!filet) {
					fprintf(stderr, "Unable to open \"%s\" for reading.\n",
							tempo);
					exit(-1);
				}

				while (!feof(filet)) {
					response = fgetc(filet);
					if (!feof(filet))
						fprintf(file2, "%c", response);
				} /* end while */

				fclose(filet);
			} /*population totale*/
		} /*boucle sur les locus*/
	} /* end for */

	if ((config.ldopt == 1) || (config.ldtotopt == 1))
		fclose(file2);

	if (config.ohtaopt == 1)
		fclose(file3);

	destroyMatald(config, *alefld);
	free(alefld);

	destroyMatald(config, *homold);
	free(homold);

	destroyMatald(config, *sumalf);
	free(sumalf);

	destroyMatald(config, *homosm);
	free(homosm);

	destroyMatdl(config, *sumld);
	free(sumld);

	destroyMatdl(config, *totsum);
	free(totsum);

	destroyMatdes(config, *burrow);
	free(burrow);

	destroyMatdes(config, *corrad);
	free(corrad);

	destroyMatst(config, *wta);
	free(wta);
}

/* ---------------Main Line Begins------------------------------------ */
int main(int argc, char **argv) {
	LinkCConfiguration config;
	int lcv;

	if (argc != 2) {
		fprintf(stderr, "USAGE:  %s <config-file>\n", argv[0]);
		return -1;
	}

	config = readConfigurationFromFile(argv[1]);
	popn = createMatpop(config);
	nmax = createMatpop(config);
	alef = createMatal(config);
	ind = createMatind(config);
	nall = createMatnumal(config);
	totn = createMatnumal(config);
	ndr1 = createMatnumal(config);
	ndr2 = createMatnumal(config);
	chisc = createMatdl(config);
	rcprob = createMatdl(config);

	/* Now we allocate the "strings" */
	aloc = (char**) malloc(sizeof(char*) * config.nlocm);
	for (lcv = 0; lcv < config.nlocm; lcv++)
		aloc[lcv] = (char*) malloc(sizeof(char) * STRING_SIZE);

	apop = (char**) malloc(sizeof(char*) * config.npopm);
	for (lcv = 0; lcv < config.npopm; lcv++)
		apop[lcv] = (char*) malloc(sizeof(char) * STRING_SIZE);

	linC(config); /* main operation occurs in here */

	destroyConfiguration(config);
	/* Clean up the memory that we dynamically allocated */
	destroyMatpop(config, popn);
	destroyMatpop(config, nmax);
	destroyMatal(config, alef);
	destroyMatind(config, ind);
	destroyMatnumal(config, nall);
	destroyMatnumal(config, totn);
	destroyMatnumal(config, ndr1);
	destroyMatnumal(config, ndr2);
	destroyMatdl(config, chisc);
	destroyMatdl(config, rcprob);

	/* Deallocate the "strings" */
	for (lcv = 0; lcv < config.nlocm; lcv++)
		free(aloc[lcv]);
	free(aloc);

	for (lcv = 0; lcv < config.npopm; lcv++)
		free(apop[lcv]);
	free(apop);

	return 0;
}
