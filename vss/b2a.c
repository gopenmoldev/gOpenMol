#include <stdio.h>

b2a(FileName)
{

      fprintf(VSS_p," 3         ! Rank value\n");
      fprintf(VSS_p," 3         ! Type of surface\n");
      fprintf(VSS_p," %d        ! number of points in the z-direction\n",npts[2]);
      fprintf(VSS_p," %d        ! number of points in the y-direction\n",npts[1]);
      fprintf(VSS_p," %d        ! number of points in the x-direction\n",npts[0]);
      fprintf(VSS_p," %f  %f    ! zmin and zmax\n",xmin[2]*autang,xmax[2]*autang);
      fprintf(VSS_p," %f  %f    ! ymin and ymax\n",xmin[1]*autang,xmax[1]*autang);
      fprintf(VSS_p," %f  %f    ! xmin and xmax\n",xmin[0]*autang,xmax[0]*autang);






