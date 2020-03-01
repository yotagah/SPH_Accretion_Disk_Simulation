#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// This program takes the a line of the data generated in the simulation and separates it accordingly to plot
int main(int argc, char **argv) {

	FILE *dados = fopen("data.dat", "r"); // Change the file name accordingly
    char *linha = (char *)malloc(1000000 * sizeof(char));
    char *split;
    int k;
    int c = 0;
    int e = atoi(argv[1]); // Number of the line wanted
    int nd = atoi(argv[2]); // Quantity of data per particle

    while(c <= e) { // Go to the wnated line
        if(fgets(linha, 1000000, dados) == NULL)
            return 1;
    	c++;
    }

    // Split the data of each particle one per line
    split = strtok(linha," ");
    while (split != NULL) {
        for(k=0; k<nd; k++) {
            printf("%g ", atof(split));
            split = strtok(NULL," ");
        }
        if(strchr(split,'\n')) {
            printf("\n");
            break;
        }
       	printf("\n");
    }

    free(linha);
    fclose(dados);

	return 0;
}
