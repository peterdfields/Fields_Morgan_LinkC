#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "LinkCUtilities.h"

#define ANSWER_LENGTH 256

char* copyString(const char *original) {
	char *ret;

	if (original == NULL )
		return copyString("");

	ret = (char*) malloc(sizeof(char) * (strlen(original) + 1));
	strcpy(ret, original);
	return ret;
}

void trimInPlace(char *string) {
	int firstNonSpace = -1;
	int lastNonSpace = -1;
	int lcv;

	for (lcv = 0; string[lcv]; lcv++) {
		if (!isspace(string[lcv])) {
			if (firstNonSpace < 0)
				firstNonSpace = lcv;

			lastNonSpace = lcv;
		}
	}

	if (firstNonSpace < 0)
		string[0] = (char) 0;
	else {
		string[lastNonSpace + 1] = (char) 0;
		memmove(string, string + firstNonSpace,
				sizeof(char) * (lastNonSpace - firstNonSpace + 2));
	}
}

int askYesNoQuestion(FILE *display, FILE *input, const char *prompt,
		int *defaultAnswer) {
	char answer[ANSWER_LENGTH];

	if (display == NULL )
		display = stdout;
	if (input == NULL )
		input = stdin;

	while (1) {
		if (prompt != NULL )
			fprintf(display, "%s ", prompt);
		if (defaultAnswer != NULL )
			fprintf(display, "(%s) ", (*defaultAnswer) ? "Y/y" : "N/n");

		if (!fgets(answer, ANSWER_LENGTH, input)) {
			fprintf(stderr, "Unable to get yes/no answer!\n");
			exit(-1);
		}

		trimInPlace(answer);
		if ((strcasecmp(answer, "yes") == 0) || (strcasecmp(answer, "y") == 0))
			return 1;
		else if ((strcasecmp(answer, "no") == 0)
				|| (strcasecmp(answer, "n") == 0))
			return 0;

		if ((answer[0] == (char) 0) && (defaultAnswer != NULL ))
			return *defaultAnswer;

		fprintf(stderr, "Please answer yes/y or no/n!\n");
		fprintf(display, "\n");
		fprintf(display, "\n");
	}

	return 0;
}
