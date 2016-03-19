#include <ctime>
#include <cmath>
#include <stdio.h>

int main(int argc, char* argv[]){
	float x, y, z;
	int N, nmr, ntypes;
	const float c = 120.27248;

	FILE *ftop;
	ftop = fopen(argv[2], "r");
	fscanf(ftop, "%5d\n", &ntypes);
	float* mass = new float [ntypes];
	int* nt = new int [ntypes];
	for(int i=0; i<ntypes; i++){
		fscanf(ftop, "%8f%5d\n", &mass[i], &nt[i]);
	}
	fclose(ftop);

	FILE *fgro;
	fgro = fopen(argv[1], "r");

	char* title = new char[5];
	fscanf(fgro, "%s\n%5d\n", title, &N);

	float* vx = new float [N];
	float* vy = new float [N];
	float* vz = new float [N];
	char tmp[5];
	for(int i=0; i<N; i++){
		fscanf(fgro, "%5d%5s%5s%5d%8f%8f%8f%8f%8f%8f\n", 
							 &nmr,tmp,tmp,&nmr,&x,&y,&z,&vx[i],&vy[i],&vz[i]);
	}
	fclose(fgro);

	float s, res = 0; int j = 0;
	for(int k=0; k<ntypes; k++){
		s = 0;
		for(int i=0; i<nt[k]; i++){
			s+=(vx[j]*vx[j]+vy[j]*vy[j]+vz[j]*vz[j]);
			j++;
		}
		s*=mass[k];
		res+=s;
	}
	res = 2*c*res/N;
	printf("Comlpete. <T> times i = %.3f\n", res);

	delete [] mass;
	delete [] nt;
	delete [] vx;
	delete [] vy;
	delete [] vz;

	return 0;
}