#ifndef __LINKDOS_CONFIG_H__
#define __LINKDOS_CONFIG_H__

typedef struct LinkCConfiguration {
	/* Data set size values */
	int npopm;
	int nlocm;
	int alm;

	/* Run configuration values */
	int ldopt;
	int ldtotopt;
	float signi;
	int ohtaopt;
	int freqopt;
	char *input;
	char *outdes;
	char *outfreq;
	char *outohta;
} LinkCConfiguration;

LinkCConfiguration readConfigurationFromFile(const char *filename);
void destroyConfiguration(LinkCConfiguration configuration);

#endif
