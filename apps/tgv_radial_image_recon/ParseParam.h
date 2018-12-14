// ParseParam.h -- See ParseParam.c for example of use

#ifndef _PARSE_PARAM
#define _PARSE_PARAM

#include <stdlib.h>
#include <iso646.h>

extern char* ReadParseParam(const char* FileName, char *VariableName);

extern char *TempPP;

#define ParseParamString(ParamFileName, Str, Len)		  \
  if ((TempPP=ReadParseParam((ParamFileName), #Str))!=NULL)	  \
    strncpy(Str, TempPP, Len); else Str[0]='\0'

#define ParseParamInt(ParamFileName, Int)			  \
  if ((TempPP=ReadParseParam((ParamFileName), #Int))!=NULL)	  \
    Int=atoi(TempPP); else Int=0

#define ParseParamHex(ParamFileName, Int)			  \
  if ((TempPP=ReadParseParam((ParamFileName), #Int))!=NULL)	  \
    Int=strtol(TempPP, NULL, 16); else Int=0

#define ParseParamFloat(ParamFileName, Flt)			  \
  if ((TempPP=ReadParseParam((ParamFileName), #Flt))!=NULL)	  \
    Flt=atof(TempPP); else Flt=0

#define ParseParamBool(ParamFileName, B)				\
  if ((TempPP=ReadParseParam((ParamFileName), #B))!=NULL)		\
    B=(toupper(TempPP[0])=='Y' || toupper(TempPP[0])=='T'|| TempPP[0]=='1'); else B=0


#endif
