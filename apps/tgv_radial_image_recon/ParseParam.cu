/******************************************************************************
  MODULE: ParseParam
  PURPOSE:  Parses a parameter file into variables
  
  EXAMPLE:
	ParseParamString("Param.txt", DirSfcBin);
	ParseParamInt("Param.txt", MaxLat);
	ParseParamFloat("Param.txt", MaxScale);
	ParseParamBool("Param.txt", DoLegend);
	ParseParamHex("Param.txt", ColorZero);
	ParseParamHex("Param.txt", ColorNan);
 	ParseParamHex("Param.txt", ColorCurve);
  will parse the following file:

-------- Start of file "Param.txt"
; This is the parameter file for RainMerge.exe and MergeMonth.exe
; Comments are preceded by ; # or !
; You can use blank lines too

DirSfcBin="/raid/data/SfcRain/"		; Use double quotes around strings
NAN=-9999.			; float, ignored value
MaxLat = 40			; integer
ColorZero=	0xC0C0C0	; hex value
ColorNan =	FFFFFF		; optional 0x
ColorCuve=	0x101040		; Syntax error
MaxScale=10.		; leading/trailing spaces are ignored
DoLegend=Y		; Bool can be T, F, Y, N, True, False, Yes, No, 0, 1...
------- End of file "Param.txt"

******************************************************************************/
#include <stdio.h>
#include <string.h>

#include "ParseParam.h"

#define LINE_DIM 1000
char *TempPP=NULL;

/******************************************************************************
  FUNCTION: ReadParam
  PURPOSE:  Read one parameter by parsing a parameter file
  RETURNS: a pointer to a string containing the value or NULL if not found
  			Use the macros to convert to typed values
******************************************************************************/
char* ReadParseParam(const char* FileName, char *VariableName) {
	static char Str[LINE_DIM];
	char *VarName, *Comment=NULL, *Equal=NULL;
	char *FirstQuote, *LastQuote, *P1, *P2;
	int Line=0, Len=0, Pos=0;
	FILE *file=fopen(FileName, "r");
	
	if (file==NULL) {
		fprintf(stderr, "\nError: Could not find file %s", FileName);
		exit(1);
	}

	while (fgets(Str, LINE_DIM-1, file) != NULL) {
		Line++;
		Len=strlen(Str);
		if (Len==0) goto Next;
		if (Str[Len-1]=='\n' or Str[Len-1]=='\r') Str[--Len]='\0';
		Equal = strchr (Str, '=');			// search for equal sign
		Pos = strcspn (Str, ";#!");			// search for comment
		Comment = (Pos==Len) ? NULL : Str+Pos;
		if (Equal==NULL or ( Comment!=NULL and Comment<=Equal)) goto Next;	// Only comment
		*Equal++ = '\0';
		if (Comment!=NULL) *Comment='\0';

		// String
		FirstQuote=strchr (Equal, '"');		// search for double quote char
		LastQuote=strrchr (Equal, '"');
		if (FirstQuote!=NULL) {
			if (LastQuote==NULL) {
				fprintf(stderr, "\nError reading parameter file %s line %d - Missing end quote.", FileName, Line);
				goto Next;
			}
			*FirstQuote=*LastQuote='\0';
			Equal=FirstQuote+1;
		}
		
		// removes leading/trailing spaces
		Pos=strspn (Str, " \t");
		if (Pos==strlen(Str)) {
			fprintf(stderr, "\nError reading parameter file %s line %d - Missing variable name.", FileName, Line);
			goto Next;		// No function name
		}
		while ((P1=strrchr(Str, ' '))!=NULL or (P2=strrchr(Str, '\t'))!=NULL)
			if (P1!=NULL) *P1='\0';
			else if (P2!=NULL) *P2='\0';
		VarName=Str+Pos;
		//while (strspn(VarName, " \t")==strlen(VarName)) VarName++;

		Pos=strspn (Equal, " \t");
		if (Pos==strlen(Equal)) {
			fprintf(stderr, "\nError reading parameter file %s line %d - Missing value.", FileName, Line);
			goto Next;		// No function name
		}
		Equal+=Pos;

//		printf("%s=%s\n", VarName, Equal);
		if (strcmp(VarName, VariableName)==0) {		// Found it
			fclose(file);
			return Equal;
		}
		Next:;
	}
	
	// not found
	fprintf(stderr, "Error reading parameter file %s - Variable %s not found.", 
				FileName, VariableName);
	fclose(file);
	return NULL;
}

