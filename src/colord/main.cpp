#include "arg_parse.h"
#include "params.h"

int main(int argc, char**argv)
{
	int params_status = parse_params(argc, argv);
	if (params_status)
		return params_status;
	
	return 0;
}