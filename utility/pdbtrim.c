/*pdbtrim  program to truncate pdb and xmol files to allow you to see what is going on in the 
center part of a block of atoms.

  The xmol part is a bit complicated because you have to keep the same atoms in each frame and in the same order, 
  and you have to write that number at the top of each frame

  Solution:  read through the first frame, enter numbers of atoms we are going to keep in keepers[i], and nkeep
  then rewind the file and start over.

  */


#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MAXI 1500
//Function Prototypes
void main(int argc, char *argv[]);


void main (int argc, char *argv[])
{
	char infilenm[40];         //temp storage for filename 
	char outfilenm[40];         //temp storage for filename 
	FILE *outfp, *infp;      	     /* pointer to output file*/
	int i, ic;
	char line[100];
	char crap[100];
	double Trunc, radius, x, y, z, dummy;
	char type, PorX;
	int natoms, ncore;
	int keepers[MAXI];
	int core[MAXI];				//atom numbers for atoms at core of viewing region
	int nkeep;
	//This next block should be enabled if you want to pick the parameter file off the command line
/*	if (argc <2) {				//no param file specified on command line
		printf("\n\n\nYOU MUST SPECIFY A PARAMETER FILE (idiot!)");
		exit(1);
	}
	else {
		while (*argv[1] != '\0'){		//load up parameter filename
			parfilenm[idum] = *argv[1];
			++idum;
			++argv[1];
		}
		outfp = readpars(); 			//read parameter file
	}  */   

//this block prompts for an input filename after execution starts
	printf ("\n\nDo you want to edit a PDB or a XMOL movie file?  [P or X]:");
	scanf ("%c", &PorX);
	PorX = toupper(PorX);


	printf ("\n\nInput File Name  BASE ONLY -- DON'T TYPE .PDB: ");
	scanf ("%s", infilenm);         //get name of parameter file
	strcpy( outfilenm, infilenm );
	switch (PorX) {
		case 'P' :
			strcat(infilenm, ".pdb");
			break;
		case 'X':
			strcat(infilenm, ".xyz");
			break;
	}


//open input file 
	if ((infp = fopen(infilenm, "r")) == NULL){
		printf("\nCannot open input file.\nMUST BE IN CURRENT DIRECTORY");
		exit (1);
	}

	printf ("\nHow many atoms are in the core of the viewing region?: ");
	scanf ("%i", &ncore);

	for (i = 0; i< ncore; ++i){
		printf ("\nAtom Number for core atom # %i: ");
		scanf ("%i", &natoms);
		core[i] = natoms;
	}

	fflush(stdin);
//Ask what kind of truncation
	printf("\nTruncate on (R)adius, (C)ylinder, or (S)quare section? [R or C or S]:");
	scanf ("%c", &type);
	type = toupper(type);

	switch( type ) {
		case 'R':
			printf("\nWhat radius (from center of face) do you want to truncate at?:");
			scanf ("%lf", &Trunc);
			switch (PorX) {
				case 'P' :
					strcat( outfilenm, "r.pdb" );		//modified files will be saved as XXXXXXr.pdb
					break;
				case 'X':
					strcat( outfilenm, "r.xyz" );		//modified files will be saved as XXXXXXr.xyz
					break;
			}
			break;
		case 'C' :
			printf("\nWhat radius do you want to truncate at?:");
			scanf ("%lf", &Trunc);
			switch (PorX) {
				case 'P' :
					strcat( outfilenm, "c.pdb" );		//modified files will be saved as XXXXXXc.pdb
					break;
				case 'X':
					strcat( outfilenm, "c.xyz" );		//modified files will be saved as XXXXXXc.xyz
					break;
			}
			break;
		case 'S' :
			printf("\nWhat edge length to truncate at?:");
			scanf ("%lf", &Trunc);
			switch (PorX) {
				case 'P' :
					strcat( outfilenm, "s.pdb" );		//modified files will be saved as XXXXXXs.pdb
					break;
				case 'X':
					strcat( outfilenm, "s.xyz" );		//modified files will be saved as XXXXXXs.xyz
					break;
			}
			break;
	}

//open the output file
	outfp = fopen(outfilenm, "w");

			/*check for succesful file open*/
	if (outfp == NULL) {
		printf("\nunable to open output file");
		exit (1);
	}



//now do it -- read file line by line, calculate radius, only write if inside the radius
	switch (PorX) {
		case 'P':
			//files are of this format:
			//  HETATM    0Mo  0         0     -15.626 -15.634  -6.188  1.00  0.00          Mo+0
			while (fgets( line, 100, infp) != NULL) { //read a line at a time !including the newline!
				for (i=0; i<100; ++i)
					crap[i] = line[i];		//copy line
//				_strnset( crap, 'a', 30 );	//overwrite the first 30 characters (makes parsing easier)
				for (i = 0; i < 30; ++i)
					crap[i] = 'a';
				sscanf (crap, "%*s%lf%lf%lf",&x, &y, &z);	//read the xyz position of atom
				switch( type ) {
					case 'R':
						radius = sqrt(x*x + y*y + z*z);
						if (radius <= Trunc)
							fputs( line, outfp );
						break;
					case 'C':
						radius = sqrt(x*x + y*y );
						if (radius <= Trunc)
							fputs( line, outfp );
						break;
					case 'S':
						if (fabs(x) <= Trunc && fabs(y) <= Trunc)
							fputs( line, outfp );
						break;
				}	//end of truncation type switch
			}		//end of while
			break;  //end of case P

		case 'X':			// a movie (xmol format) file

//first we need to read through the first frame, 
//figure out which atoms we are going to keep, enum them and put their index #s in keepers[i]
			for (i = 0; i<MAXI; ++i)
				keepers[i] = 0;
			nkeep = 0;
			fgets( line, 100, infp);		// read first line
			sscanf (line, "%i", &natoms);  //file starts with number of atoms in each frame
			fgets( line, 100, infp);		// read 2nd line
			sscanf (line, "%i", &dummy);  //line contains frame time
			for (i = 0; i < natoms; ++i){  
				fgets( line, 100, infp);		//now read atom lines
				for  (ic = 0; ic < ncore; ++ic){
					if (i == core[ic]){			// it is one of the core atoms
						++nkeep;					// so count it and add to keepers list
						keepers[i] = 1;
					}
				}
				if (keepers[i] != 1) {							//one of target atoms, so check its coordinates
					sscanf (line, "%*s%lf%lf%lf",&x, &y, &z);	//read the xyz position of atom
					switch( type ) {
						case 'R':
							radius = sqrt(x*x + y*y + z*z);
							if (radius <= Trunc){
								++nkeep;					// so count it and add to keepers list
								keepers[i] = 1;
							}
							break;
						case 'C':
							radius = sqrt(x*x + y*y );
							if (radius <= Trunc){
								++nkeep;					// so count it and add to keepers list
								keepers[i] = 1;
							}
							break;
						case 'S':
							if (fabs(x) <= Trunc && fabs(y) <= Trunc){
								++nkeep;					// so count it and add to keepers list
								keepers[i] = 1;
							}
							break;
					}	//end of truncation type switch
				}		//end of else
			}			//end of for loop

			rewind(infp);		//rewind to beginning of file
//now we are ready to process the file
			while (fgets( line, 100, infp) != NULL) { //test on the first line in each block reads, incl. the newline!
				fprintf (outfp,"%i\n", nkeep);	//write the number of atoms we are keeping
				fgets( line, 100, infp);		//next line is the frame time
				fputs (line, outfp);			// write it to the output file
				for (i = 0; i < natoms; ++i){
					fgets( line, 100, infp);		//now read atom lines
					if (keepers[i] == 1)			// on our list
						fputs (line, outfp);			// so write it to the output file
				}			//end of for loop
			}				//end of while
			break;
	}				//end of PorX switch
			
}		//end of main




