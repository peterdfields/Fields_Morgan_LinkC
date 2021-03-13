#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "LinkCConfig.h"
#include "LinkCUtilities.h"
#include "LinkCHashMap.h"

#define CONFIG_LINE_LENGTH 1024

/* The following constants are the names of the parameters from the config
 * file.
 */
#define NPOPM_KEY "npopm"
#define NLOCM_KEY "nlocm"
#define ALM_KEY "alm"
#define LDOPT_KEY "performLDComputationSubpopSubpop"
#define LDTOTOPT_KEY "performLDComputationInAllSubpopulations"
#define SIGNI_KEY "significanceLevel"
#define OHTAOPT_KEY "performOhtasLDComponentsVarianceAnalysis"
#define FREQOPT_KEY "generateIntraSubpopulationAlleleFrequenciesTable"
#define INPUT_KEY "inputDataFile"
#define OUTDES_KEY "ldOutputFile"
#define OUTFREQ_KEY "allelicFrequenciesOutputFile"
#define OUTOHTA_KEY "ohtaLDVarianceComponentsOutputFile"

static void stripComments(char *line);
static LinkCConfiguration getConfig(LinkCHashMap *map);
static const char* getRequiredParameter(LinkCHashMap *map, const char *key);
static int getBoolean(const char *key, const char *value);

LinkCConfiguration readConfigurationFromFile(const char *filename) {
	LinkCHashMap *map = createLinkCHashMap(1);
	FILE *file;
	unsigned int lineno = 0;
	char *line;
	char *value;
	int added;

	file = fopen(filename, "rt");
	if (!file) {
		fprintf(stderr, "Unable to open configuration file \"%s\".\n",
				filename);
		exit(-1);
	}

	line = (char*) malloc(sizeof(char) * CONFIG_LINE_LENGTH);
	while (fgets(line, CONFIG_LINE_LENGTH, file)) {
		added = 0;
		lineno++;
		stripComments(line);
		trimInPlace(line);
		if (line[0] == (char) 0)
			continue;
		for (value = line; !added && *value; value++) {
			if (*value == '=') {
				*value = (char) 0;
				value++;
				trimInPlace(value);
				trimInPlace(line);
				map->put(map, line, value);
				added = 1;
			}
		}

		if (!added) {
			fprintf(stderr, "Error parsing line %u of config file %s.\n",
					lineno, filename);
			exit(-1);
		}
	}
	fclose(file);

	return getConfig(map);
}

void destroyConfiguration(LinkCConfiguration configuration) {
	if (configuration.input)
		free(configuration.input);
	if (configuration.outdes)
		free(configuration.outdes);
	if (configuration.outfreq)
		free(configuration.outfreq);
	if (configuration.outohta)
		free(configuration.outohta);
}

void stripComments(char *line) {
	int lcv;

	for (lcv = 0; line[lcv]; lcv++) {
		if (line[lcv] == '#') {
			line[lcv] = (char) 0;
			break;
		}
	}
}

LinkCConfiguration getConfig(LinkCHashMap *map) {
	LinkCConfiguration config;
	memset(&config, 0, sizeof(LinkCConfiguration));

	config.npopm = strtol(getRequiredParameter(map, NPOPM_KEY), NULL, 0);
	config.nlocm = strtol(getRequiredParameter(map, NLOCM_KEY), NULL, 0);
	config.alm = strtol(getRequiredParameter(map, ALM_KEY), NULL, 0);

	config.ldopt = getBoolean(LDOPT_KEY, getRequiredParameter(map, LDOPT_KEY));
	config.ldtotopt = getBoolean(LDTOTOPT_KEY,
			getRequiredParameter(map, LDTOTOPT_KEY));

	if (config.ldopt || config.ldtotopt) {
		config.signi = strtof(getRequiredParameter(map, SIGNI_KEY), NULL );
		config.outdes = copyString(getRequiredParameter(map, OUTDES_KEY));
	}

	config.ohtaopt = getBoolean(OHTAOPT_KEY,
			getRequiredParameter(map, OHTAOPT_KEY));
	config.freqopt = getBoolean(FREQOPT_KEY,
			getRequiredParameter(map, FREQOPT_KEY));

	config.input = copyString(getRequiredParameter(map, INPUT_KEY));

	if (config.freqopt)
		config.outfreq = copyString(getRequiredParameter(map, OUTFREQ_KEY));

	if (config.ohtaopt)
		config.outohta = copyString(getRequiredParameter(map, OUTOHTA_KEY));

	return config;
}

const char* getRequiredParameter(LinkCHashMap *map, const char *key) {
	const char *ret = map->get(map, key);
	if (!ret) {
		fprintf(stderr,
				"Unable to find required configuration parameter \"%s\".\n",
				key);
		exit(-1);
	}

	return ret;
}

int getBoolean(const char *key, const char *value) {
	if (strcasecmp(value, "yes") == 0)
		return 1;
	else if (strcasecmp(value, "true") == 0)
		return 1;
	else if (strcasecmp(value, "on") == 0)
		return 1;
	else if (strcasecmp(value, "no") == 0)
		return 0;
	else if (strcasecmp(value, "false") == 0)
		return 0;
	else if (strcasecmp(value, "off") == 0)
		return 0;

	fprintf(stderr,
			"Boolean parameter %s must be yes/true/on or no/false/off!\n", key);
	exit(-1);
	return 0;
}
