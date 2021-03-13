#ifndef __LINKDOS_UTILITIES_H__
#define __LINKDOS_UTILITIES_H__

char* copyString(const char*);
void trimInPlace(char *string);
int askYesNoQuestion(FILE *display, FILE *input, const char *prompt,
		int *defaultAnswer);

#endif
