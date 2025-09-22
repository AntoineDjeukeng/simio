#include "gro.h"
#include <stdio.h>

int main(int ac, char **av)
{
	const char *json;
	FILE *fp;
	t_frame f0;
	int rc;
	if (ac < 2) {
		fprintf(stderr,"Usage: %s file.gro [channel.json]\n", av[0]);
		return 1;
	}
	if (ac >= 3) json = av[2]; else json = "data/channel.json";
	fp = fopen(av[1], "r");
	if (!fp) {
		perror("fopen");
		return 1;
	}
	rc = gro_parse(fp, &f0);
	fclose(fp);
	if (rc == -2) {
		fprintf(stderr,"Error: inconsistent atoms-per-molecule.\n");
		gro_free(&f0);
		return 1;
	}
	if (rc != 0) {
		fprintf(stderr,"Parse error.\n");
		gro_free(&f0);
		return 1;
	}
	if (gro_channel_update_or_load(&f0, json) != 0) {
		fprintf(stderr,"Note: no KIND_OTHER detected and no channel file found.\n");
	}
	gro_print_essentials(&f0);
	if (f0.sum.chan.axis != -1) {
		if (gro_channel_save(&f0, json) == 0) {
			fprintf(stderr, "Saved channel to %s\n", json);
		}
	}
	gro_free(&f0);
	return 0;
}
