#!/usr/bin/awk -f 

BEGIN {	
	ct = 1;
} 
/^>/ { 
	printf(">%s_%s\n", ctg_prefix, ct);
	++ct;
	next;
}
{
	print $0;
}