/*
		This utility evaluates the coordination numbers 
		every particle from a provided coordinates gro-file.
		
		Changes:
			added a focusing of surface particles in the input gro-file;
*/
#include <ctime>
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

#define ntype float

const ntype D_MIN = 0.5;
const ntype D_MAX = 10;
const ntype D_add = .001;

ntype min(ntype m[], int L){
	ntype tmpv (D_MAX);
	for(int i=0; i<L; i++){
		if(m[i] < tmpv) tmpv = m[i];
	}
	return tmpv;
}

int main(int argc, char* argv[]){

	int N;						// Total amoung of particles in a system;
	int devq;					// Deviation of a total charge of a system; 
	int nres, atomn, n (0);
	ntype X, Y, Z, d, r, sqr;
	ntype lx, ly, lz;			// The basic lengths of a system;
	ntype x0, y0, z0;			// Coordinates of the system center;
	ntype tmpA[] = {0, 0, 0};
	//char res[5], atom[5];

	/*argv[1] = "system.gro";
	argv[2] = "out.gro";
	argv[3] = "d_list.txt";*/
	for(int i=0; i<4;i++) 	printf("argv[%d] = %s;\n", i, argv[i]);
	if(!strcmp(argv[1],"")){
		printf(" The first argument is not defined.\n");
		return 1;
	}
	//printf("argc = %d;\n", argc);
	
	FILE *finp, *fout, *fdiam;
	if((finp = fopen(argv[1], "r"))==NULL){
		printf ("File %s is not exist!\n", argv[1]);
		return 2;
	}
	fout  = fopen(argv[2], "w");

	//Scanning of the input gro-file
	//Defining a size of a input system
	fseek(finp, -31, SEEK_END);
	fscanf(finp, "%10f%10f%10f", &lx, &ly, &lz);
	printf("A size of the system is: ( %10.5f; %10.5f; %10.5f)\n", lx, ly, lz);
	fseek(finp, 0L, SEEK_SET);

	x0 = lx/2; y0 = ly/2; z0 = lz/2;
	tmpA[0] = lx; tmpA[1] = ly; tmpA[2] = lz; 
	d  = min(tmpA, 3);				// finding a minimum from array {x0, y0, z0}
	printf("  a max available diameter is %.5f\n", d);

	char title[256];
	fgets(title, 256, finp);		//	skiping the first lines in the <finp> file
	fscanf(finp, "%5d\n", &N);
	printf("The total number of particles: N = %d\n", N);
	
	ntype* x = new ntype[N]; ntype* y = new ntype[N]; ntype* z = new ntype[N];
	ntype* vx = new ntype[N]; ntype* vy = new ntype[N]; ntype* vz = new ntype[N];
	char** res = new char* [N];
	char** atom = new char* [N];

	for(int i=0; i<N; i++){
		res[i] = new char[5];
		atom[i] = new char[5];
		fscanf(finp, "%5d%5s%5s%5d%8f%8f%8f%8f%8f%8f\n", 
			&nres,res[i],atom[i],&atomn,&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i]);
	}
	fclose(finp);
	printf("Reading the input file is complete!\n");
	
	int atn = 1;
	printf("Type a atom types number: "); scanf("%d", &atn);
	int*		 ch = new int[atn];
	char** atomt = new char*[atn];
	for(int i=0; i<atn; i++){
		atomt[i] = new char[5];
		printf(" Type name of %dth atom: ", i+1);
		scanf("%5s", atomt[i]);
		printf(" Type a charge of %s: ", atomt[i]);
		scanf("%d", &ch[i]);
	}

	if(argc > 3){
		fdiam = fopen(argv[3], "w");
		tmpA[0] = D_MAX; tmpA[1] = d; tmpA[2] = D_MAX;
		d = min( tmpA, 3); 
		fprintf(fdiam, "%s%s system file.\n   d   rescharge   atomnumbers\n",
							"The amorphous phase of TiO2 from ", argv[1]);
		printf("Start scanning of diameter values in range (%.3f, %.3f)\n", D_MIN, d);
		d = d/2; r = D_MIN/2; 
		do{
			atomn = 0; devq = 0;
			sqr = r*r;
			for(int i=0; i<N; i++){
				X = x[i]; Y = y[i]; Z = z[i];
				X -= x0; Y -= y0; Z -= z0; 
				if(X*X + Y*Y + Z*Z <= sqr){
					atomn++;
					for(int k=0; k<atn; k++){
						if(!strcmp(atom[i],atomt[k])) devq += ch[k];
					}
				}
			}
			if(devq == 0) fprintf(fdiam, "%8.3f%8d%8d\n", 2*r, devq, atomn);
			r += D_add;
		} while(r <= d);
		fclose(fdiam);
		printf("Scanning is complete! Diam list has saved in %s\n", argv[3]); 
	}

	printf("Enter a cut off diameter d: "); scanf("%f", &d);
	r = d/2; sqr = d*d/4;

	//Main function
	atomn = 0;
	fprintf(fout, "SPHERIALIZATION, d = %.2f nm, FROM %s FILE.\n%5d\n", 2*r, argv[1], 00000);
	for(int i=0; i<N; i++){
		X = x[i]; Y = y[i]; Z = z[i];
		X -= x0; Y -= y0; Z -= z0; 
		if(X*X + Y*Y + Z*Z <= sqr){
			atomn++;
			X += r; Y += r; Z += r; n++;
			fprintf(fout, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
				atomn,res[i],atom[i],atomn,x[i],y[i],z[i],vx[i],vy[i],vz[i]);
		}
		delete [] res[i]; delete atom[i];
	}
	
	fprintf(fout, "%10.5f%10.5f%10.5f\n", 2*r, 2*r, 2*r);
	fseek(fout, 0L, SEEK_SET);
	fprintf(fout, "SPHERIALIZATION, d = %.2f nm, FROM %s FILE.\n%5d\n", 2*r, argv[1], n);
	fclose(fout);

	for(int i=0; i<atn; i++) delete atomt[i];
	delete []x; delete []y; delete []z; delete []vx; delete []vy; delete []vz;
	delete [] res; delete [] atom; delete [] ch; delete [] atomt;

	printf("Complete!");

	return 0;
}