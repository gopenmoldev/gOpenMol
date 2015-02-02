/************************************************************
 * mpcube2plt.cpp
 * 
 * Molpro2002.1 cube file to gOpenMol2.1 plt file converter
 * Version 1.0.1
 *
 * Fixed bugs: 
 * - mpcube2plt couldn't read files with only one orbital or 
 *   density
 *
 * Linux/gcc:
 * g++ -lm mpcube2plt.cpp -o mpcube2plt
 *
 * Lauri Lehtovaara, lauri@st.jyu.fi
 * 8.3.2002
 ***********************************************************/

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<cstdlib>
using namespace::std;

/************************************************************
 * CONSTANTS
 ***********************************************************/
// return values
const int SUCCESS                = 0;
const int MISSING_PARAMETER      = 1;
const int INPUT_FILE_OPEN_ERROR  = 2;
const int OUTPUT_FILE_OPEN_ERROR = 3;
const int DATA_FORMAT_ERROR      = 4;
const int MEMORY_ERROR           = 5;

// data block delimiters
const string DELIMS = "\t \n";

// symbols corresponding atomic number
string ATOMSYMBOL[] = {
  "H" ,"He",
  // 2nd period
  "Li","Be","B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
  "Na","Mg","Al","Si","P" ,"S" ,"Cl","Ar",
  // 3rd period
  "K" ,"Ca",
  "Sc","Ti","V" ,"Cr","Mn","Fe","Co", "Ni","Cu","Zn",
  "Ga","Ge","As","Se","Br","Kr",
  // 4th period
  "Rb","Sr",
  "Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
  "In","Sn","Sb","Te","I" ,"Xe",
  // 5th period
  "Cs","Ba",
  "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd",
  "Tb","Dy","Ho","Er","Tm","Yb","Lu",
  "Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg",
  "Tl","Pb","Bi","Po","At","Rn",
  // 6th period
  "Fr","Ra",
  "Ac","Th","Pa","U" ,"Np","Pu","Am","Cm",
  "Bk","Cf","Es","Fm","Md","No","Lr",
  "Rf","Db","Sg","Bh","Hs","Mt"
};

double BOHR = 0.52917715; // bohr radius in angstroms


/************************************************************
 * STRUCTS
 ***********************************************************/

/************************************************************
 * cube file header
 ***********************************************************/
struct cube_header {
  string title;       // job title
  string description; // brief description of the file 
  size_t atom_count;     // number of atoms
  // origin of the grid (in bohr) 
  float grid_x_origin, grid_y_origin, grid_z_origin; 

  int n1; // 1st dimension of the grid
  // step vector for 1st dimension (in bohr)
  float n1_x_step, n1_y_step, n1_z_step; 

  int n2; // 2nd dimension of the grid
  // step vector for 1nd dimension (in bohr)
  float n2_x_step, n2_y_step, n2_z_step; 

  int n3; // 3rd dimension of the grid
  // step vector for 3rd dimension (in bohr)
  float n3_x_step, n3_y_step, n3_z_step; 

  int *atomic_number;    // list of atom numbers
  float *nuclear_charge; // list of atom charges
  // list of atom's coordinates (in bohr)
  float *atom_x_coord, *atom_y_coord, *atom_z_coord;

  bool multiple_orbs; // is there more than one orb/den 

  size_t orbital_count;   // number of orbitals
  int *orbital_number; // list of orbital numbers
};

/************************************************************
 * plt file header
 ***********************************************************/
struct plt_header {
  int rank; // rank should always be 3
  int type; // type ( 2 = orbital/density plot )
  int xdim,ydim,zdim;  // dimensions
  float xmin,ymin,zmin; // minimum coordinates
  float xmax,ymax,zmax; // maximum coordinates

  int orb_number; // number of orbital
};


/************************************************************
 * FUNCTION DECLARATIONS
 ***********************************************************/
void init_cube_header(cube_header &cubeh);
void init_plt_header(plt_header &plth);
int  read_cube_header(istream &is, cube_header &cubeh);
void print_cube_header(const cube_header &cubeh);
void print_crd_file(ostream &os, const cube_header &cubeh);
void print_plt_header(ostream &os, const cube_header &cubeh);
void write_plt_header(ostream &os, const cube_header &cubeh);
int  convert_cube_to_plt(istream &cubef,
			 string file_name_body,
			 const cube_header &cubeh);
string chop(string &,string);
string remove_whites(string);

inline int atoi(string str) { return atoi(str.c_str()); }
inline double atof(string str) { return atof(str.c_str()); }
inline string to_str(int i, int w = 0, char fc = ' ') {
  stringstream ss; 
  ss.width(w);
  ss.fill(fc);
  ss << i;
  return ss.str();
}

inline void print_col(string str, int width, ostream &os
		      , bool align_right = true) {
  os.width(width);
  if (align_right) os.setf( ios::right, ios::adjustfield ); 
  else os.setf( ios::left, ios::adjustfield );
  os << str;
}
inline void print_col(double dbl, int width, ostream &os
		      , bool align_right = true) {
  os.precision(5);
  os.width(width);
  if (align_right) os.setf( ios::right, ios::adjustfield );
  else os.setf( ios::right,ios::adjustfield );
  os.setf( ios::showpoint | ios::fixed, ios::floatfield);
  os << dbl;
}
inline void print_col(int i, int width, ostream &os,
		      bool align_right = true) {
  os.precision(5);
  os.width(width);
  if (align_right) os.setf( ios::right, ios::adjustfield ); 
  else os.setf( ios::left, ios::adjustfield );
  os << i;
}

// writes binary int to ofstream
inline void write_bin_int(ostream &os, int i_value) {
  os.write( (char*)&i_value, sizeof(i_value) );
}
// writes binary float to ofstream
inline void write_bin_float(ostream &os, float f_value) {
  os.write( (char*)&f_value, sizeof(f_value) );
}

/************************************************************
 * MAIN
 ***********************************************************/
int main(int args, char *argv[]) {
  string ofn_body; // output file name body

  if (args < 2) {
    cout << "mpcube2plt <molpro2002.1 cube file> " 
	 << "<output file name body>" << endl;
    cout << "Default output file name body is cube file "
	 << "name without file extension. " << endl
	 << "You should note that X and Z axis are swapped. "
	 << endl;
    return MISSING_PARAMETER;
  }

  if (args < 3) {
    ofn_body = argv[1];
    int pos = ofn_body.find_last_of(".");
    if ( pos != string::npos )
      ofn_body = ofn_body.substr(0,pos);
  }
  else
    ofn_body = argv[2];

  cout << "Using output file name body \""
       << ofn_body << "\"." << endl;

  // create and check input stream 
  ifstream cubefile(argv[1]);
  if (!cubefile) {
    cerr << "Couldn't read file " << argv[1] << endl;
    return INPUT_FILE_OPEN_ERROR;
  }

  // cube file header
  cube_header cubeh;
  // initialize cube header
  init_cube_header(cubeh);

  // read cube header
  int fail = read_cube_header(cubefile, cubeh);
  if ( fail == DATA_FORMAT_ERROR ) {
    cerr << "Error in cube file header." << endl;
    return fail;
  }
  else if ( fail == MEMORY_ERROR ) {
    cerr << "Failed to allocate memory." << endl;
    return fail;
  }

  cout << "----------------------------------------" << endl;
  cout << "Cube header:" << endl;
  cout << "----------------------------------------" << endl;
  print_cube_header(cubeh);

  cout << "----------------------------------------" << endl;
  cout << "Crd file:" << endl;
  cout << "----------------------------------------" << endl;
  print_crd_file(cout, cubeh);

  ofstream crdf( (ofn_body + ".crd").c_str() );
  print_crd_file(crdf, cubeh);
  crdf.close();

  cout << "----------------------------------------" << endl;
  cout << "Plt file header:" << endl;
  cout << "----------------------------------------" << endl;
  print_plt_header(cout, cubeh);

  cout << "----------------------------------------" << endl;
  cout << "Converting cube file to plt files..." << endl;
  cout << "----------------------------------------" << endl;

  fail = convert_cube_to_plt(cubefile, ofn_body, cubeh);
  if ( fail == OUTPUT_FILE_OPEN_ERROR ) {
    cerr << endl << "Couldn't create output file." << endl;
    return fail;
  }
  else if ( fail == DATA_FORMAT_ERROR ) {
    cerr << endl << "Error in cube file data." << endl;
  }

  cout << "----------------------------------------" << endl;
  cout << "Conversion complete..." << endl;
  cout << "----------------------------------------" << endl;

  return SUCCESS;
}

/************************************************************
 * Initializes cube file header
 ***********************************************************/
void init_cube_header(cube_header &cubeh) {
  cubeh.atomic_number  = (int*) 0;
  cubeh.nuclear_charge = (float*) 0;
  cubeh.atom_x_coord   = (float*) 0;
  cubeh.atom_y_coord   = (float*) 0;
  cubeh.atom_z_coord   = (float*) 0;
  cubeh.orbital_number = (int*) 0;
}

/************************************************************
 * initialize plt file header
 ***********************************************************/
void init_plt_header(plt_header &plth) {
  plth.rank = 3; // should always be 3
  plth.type = 2; // orbital/density plot
}

/************************************************************
 * Reads and checks cube file header
 ***********************************************************/
int read_cube_header(istream &is, cube_header &cubeh) {
  size_t i;
  string str;

  // read 1st line
  if ( ! getline(is,str) ) 
    return DATA_FORMAT_ERROR; 
  // chop title
  cubeh.title = chop(str,DELIMS);

  // read 2nd line
  if ( ! getline(is,str) ) 
    return DATA_FORMAT_ERROR; 
  // chop description
  cubeh.description = chop(str,DELIMS);

  // read 3st line
  if ( ! getline(is,str) ) 
    return DATA_FORMAT_ERROR; 
  // chop number of atoms
  int atom_count = atoi( chop(str,DELIMS) );
  cubeh.atom_count = abs( atom_count );
  // if atom count is negative, then there is more 
  // than one orbital
  if ( atom_count < 0 )
    cubeh.multiple_orbs = true;
  else
    cubeh.multiple_orbs = false;

  // chop grid origin
  cubeh.grid_x_origin = atof( chop(str,DELIMS) );
  cubeh.grid_y_origin = atof( chop(str,DELIMS) );
  cubeh.grid_z_origin = atof( chop(str,DELIMS) );

  // read 4th line
  if ( ! getline(is,str) ) 
    return DATA_FORMAT_ERROR; 
  // chop number grid points in 1st dimension
  cubeh.n1 = abs( atoi( chop(str,DELIMS) ) );
  // chop n1 step
  cubeh.n1_x_step = atof( chop(str,DELIMS) );
  cubeh.n1_y_step = atof( chop(str,DELIMS) );
  cubeh.n1_z_step = atof( chop(str,DELIMS) );

  // read 5th line
  if ( ! getline(is,str) ) 
    return DATA_FORMAT_ERROR; 
  // chop number grid points in 2nd dimension
  cubeh.n2 = abs( atoi( chop(str,DELIMS) ) );
  // chop n2 step
  cubeh.n2_x_step = atof( chop(str,DELIMS) );
  cubeh.n2_y_step = atof( chop(str,DELIMS) );
  cubeh.n2_z_step = atof( chop(str,DELIMS) );

  // read 6th line
  if ( ! getline(is,str) ) 
    return DATA_FORMAT_ERROR; 
  // chop number grid points in 3rd dimension
  cubeh.n3 = abs( atoi( chop(str,DELIMS) ) );
  // chop n3 step
  cubeh.n3_x_step = atof( chop(str,DELIMS) );
  cubeh.n3_y_step = atof( chop(str,DELIMS) );
  cubeh.n3_z_step = atof( chop(str,DELIMS) );

  cout << "test" << cubeh.atom_count <<endl;

  // do we have atoms?
  if ( cubeh.atom_count > 0 ) {
    // if we do, allocate some memory
    cubeh.atomic_number  = new int[cubeh.atom_count];
    cubeh.nuclear_charge = new float[cubeh.atom_count];
    cubeh.atom_x_coord   = new float[cubeh.atom_count];
    cubeh.atom_y_coord   = new float[cubeh.atom_count];
    cubeh.atom_z_coord   = new float[cubeh.atom_count];
    
    // and check if success
    if ( ( ! cubeh.atomic_number  ) ||
	 ( ! cubeh.nuclear_charge ) ||
	 ( ! cubeh.atom_x_coord   ) ||
	 ( ! cubeh.atom_y_coord   ) ||
	 ( ! cubeh.atom_z_coord   ) ) {
      return MEMORY_ERROR;
    }      
  }

  // read atomic data
  for(i = 0; i < cubeh.atom_count; i++) {
    // read (7+i)th line
    if ( ! getline(is,str) ) 
      return DATA_FORMAT_ERROR; 
    // chop atomic numer
    cubeh.atomic_number[i] = 
      abs( atoi( chop(str,DELIMS) ) );
    // chop nuclear charge
    cubeh.nuclear_charge[i] = atof( chop(str,DELIMS) );
    // chop atom's coordinates
    cubeh.atom_x_coord[i] = atof( chop(str,DELIMS) );
    cubeh.atom_y_coord[i] = atof( chop(str,DELIMS) );
    cubeh.atom_z_coord[i] = atof( chop(str,DELIMS) );
  }

  // if more than one orbital or density
  if ( cubeh.multiple_orbs ) {
    
    // read (8+atom_count)th line
    if ( ! getline(is,str) ) 
      return DATA_FORMAT_ERROR; 
    // chop number of orbitals
    cubeh.orbital_count = abs( atoi( chop(str,DELIMS) ) );
  
  
    // do we have orbitals or densities?
    if (cubeh.orbital_count > 0) {
      // if do, allocate memory
      cubeh.orbital_number = new int[cubeh.orbital_count];
    
      // check if success
      if ( ! cubeh.orbital_number ) {
	return MEMORY_ERROR;
      }
    }
  
      // chop orbital numbers
      for(i = 0; i < cubeh.orbital_count; i++) {
	// if no more orbital numbers on this line
	if (str == "")
	  // read line
	  if ( ! getline(is,str) ) 
	    return DATA_FORMAT_ERROR; 
      
	// chop orbital number
	cubeh.orbital_number[i] =
	  abs( atoi( chop(str,DELIMS) ) );
      }
  }
  else {
    cubeh.orbital_count = 1;
    cubeh.orbital_number = new int[1];
    cubeh.orbital_number[0] = 1;
  }

  return SUCCESS;
}

/************************************************************
 * Prints cube file header
 ***********************************************************/
void print_cube_header(const cube_header &cubeh) {
  size_t i;
  
  cout << "Number of atoms: " << cubeh.atom_count << endl;
  cout << "Grid origin: " 
       << cubeh.grid_x_origin << ", "
       << cubeh.grid_x_origin << ", "
       << cubeh.grid_x_origin << endl;

  cout << "Dimensions: " << cubeh.n1 << ", "
       << cubeh.n2 << ", " << cubeh.n3 << endl;
  cout << "N1 step vector: "
       << cubeh.n1_x_step << ", "
       << cubeh.n1_y_step << ", "
       << cubeh.n1_z_step << endl;
  cout << "N2 step vector: "
       << cubeh.n2_x_step << ", "
       << cubeh.n2_y_step << ", "
       << cubeh.n2_z_step << endl;
  cout << "N3 step vector: "
       << cubeh.n3_x_step << ", "
       << cubeh.n3_y_step << ", "
       << cubeh.n3_z_step << endl;
  
  for(i = 0; i < cubeh.atom_count; i++) {
    cout << "Atom " << (i+1) << ":" << endl;
    cout << "\tAtomic number:  " 
	 << cubeh.atomic_number[i] << endl;
    cout << "\tNuclear charge: " 
	 << cubeh.nuclear_charge[i] << endl;
    cout << "\tAtomic coordinates: " 
	 << cubeh.atom_x_coord[i] << ", "
	 << cubeh.atom_y_coord[i] << ", "
	 << cubeh.atom_z_coord[i] << endl;
  }

  for(i = 0; i < cubeh.orbital_count; i++)
    cout << "Orbital " << (i+1) << " orbital number: " 
	 << cubeh.orbital_number[i] << endl;
}

/************************************************************
 * Prints crd file to stream from cube file header
 ***********************************************************/
void print_crd_file(ostream &os, const cube_header &cubeh) {
  size_t i;

  os << "* " << cubeh.title << endl;

  // number of atoms
  print_col((int)cubeh.atom_count,5,os);
  os << endl;
  
  for(i = 0; i < cubeh.atom_count; i++) {
    // atom id 
    print_col((int)i+1,5,os);
    os<<' ';
    // resno
    print_col("1",4,os);
    os<<' ';
    // restype
    print_col("GAUS",4,os);
    os<<' ';
    // atom symbol
    print_col( ATOMSYMBOL[ cubeh.atomic_number[i] - 1 ], 
	       4, os, false);
    
    //------------------------------
    // NOTE: X AND Z SWAPPED
    //------------------------------
    // x coordinate
    print_col( cubeh.atom_z_coord[i] * BOHR, 10, os);
    // y coordinate
    print_col( cubeh.atom_y_coord[i] * BOHR, 10, os);
    // z coordinate
    print_col( cubeh.atom_x_coord[i] * BOHR, 10, os);
    
    os<<' ';
    // segid
    print_col("GAUS",4,os);
    os<<' ';
    // resid
    print_col("1",4,os,false);

    print_col(0.0,10,os);

    os << endl;
  }

}

/************************************************************
 * Prints plt file header to stream from cube file header
 ***********************************************************/
void print_plt_header(ostream &os, 
		      const cube_header &cubeh) {

  // print rank (3 always) and type 2 (orbital/density plot)
  os << "3 2" << endl;
  // dimensions
  os << cubeh.n1 << ' ' << cubeh.n2 << ' ' << cubeh.n3 
     << endl;

  // bounds of coordinates
  //------------------------------
  // NOTE: X AND Z SWAPPED
  //------------------------------

  // zmin zmax
  os << cubeh.grid_x_origin * BOHR << ' '
     << - cubeh.grid_x_origin * BOHR << endl;
  // ymin ymax
  os << cubeh.grid_y_origin * BOHR << ' '
     << - cubeh.grid_y_origin * BOHR << endl;
  // xmin xmax
  os << cubeh.grid_z_origin * BOHR << ' '
     << - cubeh.grid_z_origin * BOHR << endl;
}

/************************************************************
 * Writes plt file header to binary file from cube file 
 * header
 ***********************************************************/
void write_plt_header(ostream &os, 
		      const cube_header &cubeh) {
  // print rank (3 always) and type 2 (orbital/density plot)
  write_bin_int( os, 3 );
  write_bin_int( os, 2 );
  // dimensions
  write_bin_int( os, cubeh.n1 );
  write_bin_int( os, cubeh.n2 ); 
  write_bin_int( os, cubeh.n3 );

  // bounds of coordinates
  //------------------------------
  // NOTE: X AND Z SWAPPED
  //------------------------------

  // zmin zmax
  write_bin_float( os, (float)( cubeh.grid_x_origin * BOHR ) );
  write_bin_float( os, (float)( - cubeh.grid_x_origin * BOHR ) );
  // ymin ymax
  write_bin_float( os, (float)( cubeh.grid_y_origin * BOHR ) );
  write_bin_float( os, (float)( - cubeh.grid_y_origin * BOHR ) );
  // xmin xmax
  write_bin_float( os, (float)( cubeh.grid_z_origin * BOHR ) );
  write_bin_float( os, (float)( - cubeh.grid_z_origin * BOHR ) );
}

/************************************************************
 * Converts cube file data to multiple plt files
 ***********************************************************/
int convert_cube_to_plt(istream &cubef, 
			string file_name_body,
			const cube_header &cubeh) {

  int    i = 0, j = 0, k = 0;
  int    c = 0;      // counter when need to read more 
  size_t n = 0;
  string str = "";

  /* molpro format:
   * Loops from outtermost to innermost: n1,n2,orbs,n3.
   * So there is n1*n2 blocks, which contains orbital's
   * density information. In each block there is number
   * of orbs subblocks, which are n3 long. The end of
   * the block (not subblocks) are filled with zeros
   * (if needed), so that blocks are multiples of six.
   */

  // initialize and print header to plt files
  cout << "Creating output files...";

  string file_name;

  ofstream *pltf = new ofstream[cubeh.orbital_count];

  for(n = 0; n < cubeh.orbital_count; n++) {

    // opening file
    file_name = file_name_body 
      + to_str( cubeh.orbital_number[n], 3, '_');
    pltf[n].open( (file_name + ".plt").c_str(), 
		  ios::binary | ios::trunc | ios::out );

    // checking if success
    if ( ! pltf[n].is_open() )
      return OUTPUT_FILE_OPEN_ERROR;

    if ( n % 5 == 0 ) cout << endl; 
    cout << file_name << ".plt  ";

    write_plt_header(pltf[n],cubeh);
  }
  cout << endl;

  cout << "Parsing blocks (" 
       << cubeh.n1*cubeh.n2 << " total): ";

  // loop n1
  for(i = 0; i < cubeh.n1; i++) {

    if (i % 10 == 0) 
      cout << endl;
    cout.width(5);
    cout << i*cubeh.n2 << "  ";

    // loop n2
    for(j = 0; j < cubeh.n2; j++) {
      // discard remainings of previous data block
      c   = 0;
      str = "";

      // loop orbs
      for(n = 0; n < cubeh.orbital_count; n++) {
	// loop n3
	for(k = 0; k < cubeh.n3; k++) {
      
	  // if six is chopped, read new data line
	  if (c % 6 == 0 ) {
	    if ( ! getline(cubef,str) )
	      return DATA_FORMAT_ERROR;
	    // prevent overflow
	    c = 0;
	  }
	  // increase counter read counter
	  c++;

	  // write data to orbital files 
	  write_bin_float( pltf[n], (float)( atof( chop(str,DELIMS) ) ) );
	}
      }
    }
  }

  // close files and free resources
  for(n = 0; n < cubeh.orbital_count; n++)
    pltf[n].close();

  delete[] pltf;

  cout << endl;

  return SUCCESS;
}



/************************************************************
 * strmodif.cpp file included here
 ***********************************************************/

bool is_white(char c) {
       if ( c == '\t' ) return true;
  else if ( c == ' '  ) return true;
  else if ( c == '\n'  ) return true;
  else return false;
}

string remove_whites(string str) {
  string::size_type i;
  for(i = 0; is_white(str[i]) && i < str.size(); i++ );
  return str.substr(i,string::npos);
}

/************************************************************
 * Chops from wood until delimeter then returns chopped logs
 ***********************************************************/
string chop(string &wood, string delims) {
  string log;
  char c;
  string::size_type i;

  wood = remove_whites(wood);
  i = wood.find_first_of(delims);

  if (i == string::npos) {
    log = wood;
    wood = "";
  }
  else {
    log = wood.substr(0,i);
    wood = wood.substr(i+1,string::npos);
  }

  return log;
}

