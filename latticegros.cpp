/*
 It takes two arguments: both paths of a input file and 
 a output file.
 */
#include <ctime>
#include <cmath>
#include <stdio.h>

#define ntype float

int main(int argc, char* argv[]){
//	int c = 0; 
	ntype x, y, z, x0, y0, z0,
			xn, yn, zn, vx, vy, vz, d;
			vx = 0; vy = 0; vz = 0;
	int na, nb, nc;

	FILE *finp, *fgro;
	finp = fopen(argv[1], "r");

	ntype mvec[3] = { 0, 0, 0};	// the main vectors of a unit cell
	int atn = 0;						// a atom types number 
	char title[256];					// a title of a input file
	char tmp[5];

	//Reading a unit cell configuration file 
	fgets(title, 256, finp);
	fscanf(finp, "%8f%8f%8f%5s\n", &mvec[0], &mvec[1], &mvec[2], tmp);
	fscanf(finp, "%3d\n", &atn);
	mvec[0] = mvec[0]/10; mvec[1] = mvec[1]/10; mvec[2] = mvec[2]/10; 

	char** atom = new char*[atn];	// names of atom types 
	char** res  = new char*[atn];	// a residue name of the whole of a crystal lattice
	int*	tpu = new int[atn];		// the atom number of every type per unit	

	int n = 0;
	for(int i=0; i<atn; i++){
		atom[i] = new char[5];
		fscanf(finp, "%5s%5d", atom[i], &tpu[i]);
		n += tpu[i];
	}
	int ppu = n;					// a number of particles per unit
	
	ntype*** coord = new ntype**[atn];
	
	//int l = -1;
	for(int i=0; i<atn; i++){
		coord[i] = new ntype*[tpu[i]];
		res[i]  = new char[5];
		fscanf(finp, "\n%5s", res[i]);
		for(int k=0; k<tpu[i]; k++){
			//l++;
			coord[i][k] = new ntype[3];
			fscanf(finp, "\n%8f%8f%8f", &coord[i][k][0], &coord[i][k][1], &coord[i][k][2]);
		}
	}
	fclose(finp);
	
	//printf("Reading of input file is complete!\nType the diameter, nm: ");
	scanf("%f", &d);
	ntype sqr = d*d/4;		// square of raduis
	na = ceil(d/mvec[0]);
	nb = ceil(d/mvec[1]);
	nc = ceil(d/mvec[2]);

	//Forming of the output .gro file
	fgro = fopen(argv[2], "w");
	int N = ppu*na*nb*nc; int t = time(0);
	fprintf(fgro, "SPHERICAL PARTICLE OF TiO2 (%s), D = %.3f nm, %dx%dx%d in %5d\n%5d\n", argv[3], d, na, nb, nc, t, 00000);
	
	int natom = 0;
	int nres = 0;
	x0 = mvec[0]*na/2; y0 = mvec[1]*nb/2; z0 = mvec[2]*nc/2;
	for(int ktype=0; ktype<atn; ktype++){
		for(int k=0; k<nc; k++){
			zn = k*mvec[2];
			for(int j=0; j<nb; j++){
				yn = j*mvec[1];
				for(int i=0; i<na; i++){
					xn = i*mvec[0];
					for(int k=0; k<tpu[ktype]; k++){
						x = xn + coord[ktype][k][0]/10 - x0;
						y = yn + coord[ktype][k][1]/10 - y0;
						z = zn + coord[ktype][k][2]/10 - z0;
						if (x*x + y*y + z*z <= sqr){
							x += x0; y += y0; z += z0; natom++;
							fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", 
										natom,res[ktype],atom[ktype],natom,x,y,z,vx,vy,vz);
						}
					}
				}
			}
		}
	}
	//na++; nb++; nc++;
	fprintf(fgro, "%10.5f%10.5f%10.5f\n", mvec[0]*na, mvec[1]*nb, mvec[2]*nc);
	fseek(fgro, 0L, SEEK_SET);
	fprintf(fgro, "SPHERICAL PARTICLE OF TiO2 (%s), D = %.3f nm, %dx%dx%d in %5d\n%5d", 
			  argv[3], d, na, nb, nc, t, natom);	
	fclose(fgro);
	//printf("Comlpete. Amount = %d atoms\n", natom);
	printf("%d\n", natom);

	for(int i=0; i<atn; i++){
		delete [] atom[i];
		delete [] res[i];
		for(int k=0; k<tpu[i]; k++){
			delete [] coord[i][k];
		}
		delete [] coord[i];
	}
	delete [] tpu; delete [] atom; delete [] res; delete [] coord;

	return 0;
}