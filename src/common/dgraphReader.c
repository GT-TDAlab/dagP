/*
 *   From Matrix Market I/O library for ANSI C
 *
 *   See http://math.nist.gov/MatrixMarket for details.
 *
 *
 */
/*
 ***************************************************************
 *   modifications by Bora Ucar, and others.                   *
 ***************************************************************
 *  symmetric case is handled in such a way that, I and J      *
 *            now contains both symmetric entries, nz is       *
 *            updated accordingly.                             *
 ***************************************************************
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <fcntl.h>

#include "dgraphReader.h"

#define PRINT_PROGRESS


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wunused-parameter"


#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64


typedef char MM_typecode[4];

char *mm_typecode_to_str(MM_typecode matcode);


/********************* MM_typecode query fucntions ***************************/

#define mm_is_matrix(typecode)  ((typecode)[0]=='M')

#define mm_is_sparse(typecode)  ((typecode)[1]=='C')
#define mm_is_coordinate(typecode)((typecode)[1]=='C')
#define mm_is_dense(typecode)   ((typecode)[1]=='A')
#define mm_is_array(typecode)   ((typecode)[1]=='A')

#define mm_is_complex(typecode) ((typecode)[2]=='C')
#define mm_is_real(typecode)        ((typecode)[2]=='R')
#define mm_is_pattern(typecode) ((typecode)[2]=='P')
#define mm_is_integer(typecode) ((typecode)[2]=='I')

#define mm_is_symmetric(typecode)((typecode)[3]=='S')
#define mm_is_general(typecode) ((typecode)[3]=='G')
#define mm_is_skew(typecode)    ((typecode)[3]=='K')
#define mm_is_hermitian(typecode)((typecode)[3]=='H')

int mm_is_valid(MM_typecode matcode);       /* too complex for a macro */


/********************* MM_typecode modify fucntions ***************************/

#define mm_set_matrix(typecode) ((*typecode)[0]='M')
#define mm_set_coordinate(typecode) ((*typecode)[1]='C')
#define mm_set_array(typecode)  ((*typecode)[1]='A')
#define mm_set_dense(typecode)  mm_set_array(typecode)
#define mm_set_sparse(typecode) mm_set_coordinate(typecode)

#define mm_set_complex(typecode)((*typecode)[2]='C')
#define mm_set_real(typecode)   ((*typecode)[2]='R')
#define mm_set_pattern(typecode)((*typecode)[2]='P')
#define mm_set_integer(typecode)((*typecode)[2]='I')


#define mm_set_symmetric(typecode)((*typecode)[3]='S')
#define mm_set_general(typecode)((*typecode)[3]='G')
#define mm_set_skew(typecode)   ((*typecode)[3]='K')
#define mm_set_hermitian(typecode)((*typecode)[3]='H')

#define mm_clear_typecode(typecode) ((*typecode)[0]=(*typecode)[1]= \
(*typecode)[2]=' ',(*typecode)[3]='G')

#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)


/********************* Matrix Market error codes ***************************/


#define MM_COULD_NOT_READ_FILE  11
#define MM_PREMATURE_EOF        12
#define MM_NOT_MTX              13
#define MM_NO_HEADER            14
#define MM_UNSUPPORTED_TYPE     15
#define MM_LINE_TOO_LONG        16
#define MM_COULD_NOT_WRITE_FILE 17


/******************** Matrix Market internal definitions ********************

 MM_matrix_typecode: 4-character sequence

 ojbect         sparse/     data        storage
 dense      type        scheme

 string position:    [0]        [1]         [2]         [3]

 Matrix typecode:  M(atrix)  C(oord)        R(eal)      G(eneral)
 A(array)   C(omplex)   H(ermitian)
 P(attern)   S(ymmetric)
 I(nteger)  K(kew)

 ***********************************************************************/

#define MM_MTX_STR      "matrix"
#define MM_ARRAY_STR    "array"
#define MM_DENSE_STR    "array"
#define MM_COORDINATE_STR "coordinate"
#define MM_SPARSE_STR   "coordinate"
#define MM_COMPLEX_STR  "complex"
#define MM_REAL_STR     "real"
#define MM_INT_STR      "integer"
#define MM_GENERAL_STR  "general"
#define MM_SYMM_STR     "symmetric"
#define MM_HERM_STR     "hermitian"
#define MM_SKEW_STR     "skew-symmetric"
#define MM_PATTERN_STR  "pattern"



int mm_is_valid(MM_typecode matcode)
{
    if (!mm_is_matrix(matcode)) return 0;
    if (mm_is_dense(matcode) && mm_is_pattern(matcode)) return 0;
    if (mm_is_real(matcode) && mm_is_hermitian(matcode)) return 0;
    if (mm_is_pattern(matcode) && (mm_is_hermitian(matcode) ||
                                   mm_is_skew(matcode))) return 0;
    return 1;
}

int mm_read_banner(FILE *f, MM_typecode *matcode)
{
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH];
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;
    /*    int ret_code;*/


    mm_clear_typecode(matcode);

    if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
        return MM_PREMATURE_EOF;

    if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type,
               storage_scheme) != 5)
        return MM_PREMATURE_EOF;

    for (p=mtx; *p!='\0'; *p = (char) tolower(*p),p++);  /* convert to lower case */
    for (p=crd; *p!='\0'; *p = (char) tolower(*p),p++);
    for (p=data_type; *p!='\0'; *p = (char) tolower(*p),p++);
    for (p=storage_scheme; *p!='\0'; *p = (char) tolower(*p),p++);

    /* check for banner */
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
        return MM_NO_HEADER;

    /* first field should be "mtx" */
    if (strcmp(mtx, MM_MTX_STR) != 0)
        return  MM_UNSUPPORTED_TYPE;
    mm_set_matrix(matcode);


    /* second field describes whether this is a sparse matrix (in coordinate
     storgae) or a dense array */


    if (strcmp(crd, MM_SPARSE_STR) == 0)
        mm_set_sparse(matcode);
    else
    if (strcmp(crd, MM_DENSE_STR) == 0)
        mm_set_dense(matcode);
    else
        return MM_UNSUPPORTED_TYPE;


    /* third field */

    if (strcmp(data_type, MM_REAL_STR) == 0)
        mm_set_real(matcode);
    else
    if (strcmp(data_type, MM_COMPLEX_STR) == 0)
        mm_set_complex(matcode);
    else
    if (strcmp(data_type, MM_PATTERN_STR) == 0)
        mm_set_pattern(matcode);
    else
    if (strcmp(data_type, MM_INT_STR) == 0)
        mm_set_integer(matcode);
    else
        return MM_UNSUPPORTED_TYPE;


    /* fourth field */

    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
        mm_set_general(matcode);
    else
    if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
        mm_set_symmetric(matcode);
    else
    if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        mm_set_hermitian(matcode);
    else
    if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
        mm_set_skew(matcode);
    else
        return MM_UNSUPPORTED_TYPE;


    return 0;
}

int mm_write_mtx_crd_size(FILE *f, idxType M, idxType N, idxType nz)
{
    if(sizeof(idxType) == sizeof(int))
    {
        if (fprintf(f, "%d %d %d\n", M, N, nz) != 3)
            return MM_COULD_NOT_WRITE_FILE;
        else
            return 0;
    }
    else
    {
        if (fprintf(f, "%ld %ld %ld\n", M, N, nz) != 3)
            return MM_COULD_NOT_WRITE_FILE;
        else
            return 0;
    }
}

int mm_read_mtx_crd_size(FILE *f, idxType *M, idxType *N, idxType *nz )
{
    char line[MM_MAX_LINE_LENGTH];
    /*    int ret_code;*/
    int num_items_read;

    /* set return null parameter values, in case we exit with errors */
    *M = *N = *nz = 0;

    /* now continue scanning until you reach the end-of-comments */
    do
    {
        if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL)
            return MM_PREMATURE_EOF;
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if (sizeof(idxType) == sizeof(int))
    {
        if (sscanf(line, "%d %d %d", M, N, nz) == 3)
            return 0;

        else
            do
            {
                num_items_read = fscanf(f, "%d %d %d", M, N, nz);
                if (num_items_read == EOF) return MM_PREMATURE_EOF;
            }
            while (num_items_read != 3);
    }
    else
    {
        if (sscanf(line, "%ld %ld %ld", M, N, nz) == 3)
            return 0;

        else
            do
            {
                num_items_read = fscanf(f, "%ld %ld %ld", M, N, nz);
                if (num_items_read == EOF) return MM_PREMATURE_EOF;
            }
            while (num_items_read != 3);
    }
    return 0;
}


int mm_read_mtx_array_size(FILE *f, idxType *M, idxType *N)
{
    char line[MM_MAX_LINE_LENGTH];
    /*    int ret_code;*/
    int num_items_read;

    /* set return null parameter values, in case we exit with errors */
    *M = *N = 0;

    /* now continue scanning until you reach the end-of-comments */
    do
    {
        if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL)
            return MM_PREMATURE_EOF;
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if(sizeof(idxType) == sizeof(int))
    {
        if (sscanf(line, "%d %d", M, N) == 2)
            return 0;

        else /* we have a blank line */
            do
            {
                num_items_read = fscanf(f, "%d %d", M, N);
                if (num_items_read == EOF) return MM_PREMATURE_EOF;
            }
            while (num_items_read != 2);
    }
    else
    {
        if (sscanf(line, "%ld %ld", M, N) == 2)
            return 0;

        else /* we have a blank line */
            do
            {
                num_items_read = fscanf(f, "%ld %ld", M, N);
                if (num_items_read == EOF) return MM_PREMATURE_EOF;
            }
            while (num_items_read != 2);

    }
    return 0;
}

int mm_write_mtx_array_size(FILE *f, idxType M, idxType N)
{
    if(sizeof(idxType) == sizeof(int))
    {
        if (fprintf(f, "%d %d\n", M, N) != 2)
            return MM_COULD_NOT_WRITE_FILE;
        else
            return 0;
    }
    else
    {
        if (fprintf(f, "%ld %ld\n", M, N) != 2)
            return MM_COULD_NOT_WRITE_FILE;
        else
            return 0;
    }
}



/*-------------------------------------------------------------------------*/

/******************************************************************/
/* use when I[], J[], and val[]J, and val[] are already allocated */
/******************************************************************/

int mm_read_mtx_crd_data(FILE *f, idxType M, idxType N, idxType nz, idxType I[], idxType J[],
                         double val[], MM_typecode matcode)
{
    idxType i;
    if (mm_is_complex(matcode))
    {
        for (i=0; i<nz; i++)
        {
            if(sizeof(idxType) == sizeof(int))
            {
                if (fscanf(f, "%d %d %lg %lg", &I[i], &J[i], &val[2*i], &val[2*i+1])
                    != 4) return MM_PREMATURE_EOF;
            }
            else
            {
                if (fscanf(f, "%ld %ld %lg %lg", &I[i], &J[i], &val[2*i], &val[2*i+1])
                    != 4) return MM_PREMATURE_EOF;
            }
        }
    }
    else if (mm_is_real(matcode))
    {
        for (i=0; i<nz; i++)
        {
            if(sizeof(idxType) == sizeof(int))
            {
                if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i])
                    != 3) return MM_PREMATURE_EOF;
            }
            else
            {
                if (fscanf(f, "%ld %ld %lg\n", &I[i], &J[i], &val[i])
                    != 3) return MM_PREMATURE_EOF;
            }

        }
    }

    else if (mm_is_pattern(matcode))
    {
        for (i=0; i<nz; i++)
        {
            if (sizeof(idxType) == sizeof(int)){
                if (fscanf(f, "%d %d", &I[i], &J[i])!= 2) return MM_PREMATURE_EOF;
            }
            else
            if (fscanf(f, "%ld %ld", &I[i], &J[i])!= 2) return MM_PREMATURE_EOF;

        }
    }
    else if (mm_is_integer(matcode))
    {
        for (i=0; i<nz; i++)
        {
            if(sizeof(idxType) == sizeof(int))
            {
                if (fscanf(f, "%d %d %lf\n", &I[i], &J[i], &val[i])
                    != 3) return MM_PREMATURE_EOF;
            }
            else
            {
                if (fscanf(f, "%ld %ld %lf\n", &I[i], &J[i], &val[i])
                    != 3) return MM_PREMATURE_EOF;
            }
        }
    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;

}

int mm_read_mtx_crd_entry(FILE *f, idxType *I, idxType *J,
                          double *real, double *imag, MM_typecode matcode)
{
    /*    int i;*/
    if (mm_is_complex(matcode))
    {
        if(sizeof(I[0]) == sizeof(int))
        {
            if (fscanf(f, "%d %d %lg %lg", I, J, real, imag)
                != 4) return MM_PREMATURE_EOF;
        }
        else
        {
            if (fscanf(f, "%ld %ld %lg %lg", I, J, real, imag)
                != 4) return MM_PREMATURE_EOF;
        }
    }
    else if (mm_is_real(matcode))
    {
        if(sizeof(I[0]) == sizeof(int))
        {
            if (fscanf(f, "%d %d %lg\n", I, J, real)
                != 3) return MM_PREMATURE_EOF;
        }
        else
        if (fscanf(f, "%ld %ld %lg\n", I, J, real)
            != 3) return MM_PREMATURE_EOF;
    }

    else if (mm_is_pattern(matcode))
    {
        if(sizeof(I[0]) == sizeof(int))
        {
            if (fscanf(f, "%d %d", I, J) != 2) return MM_PREMATURE_EOF;
        }
        else
        if (fscanf(f, "%ld %ld", I, J) != 2) return MM_PREMATURE_EOF;

    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;

}

int mm_num_all_nnz(char *fname)
{
    int ret_code;
    idxType M, N, nz;
    idxType add = 0;
    MM_typecode matcode;

    FILE *f;
    if (strcmp(fname, "stdin") == 0) f=stdin;
    else
    if ((f = fopen(fname, "r")) == NULL)
        return MM_COULD_NOT_READ_FILE;


    if ((ret_code = mm_read_banner(f, &matcode)) != 0)
        return ret_code;

    if (!(mm_is_valid(matcode) && mm_is_sparse(matcode) &&
          mm_is_matrix(matcode)))
        return MM_UNSUPPORTED_TYPE;
    if(mm_is_complex(matcode))
        return MM_UNSUPPORTED_TYPE;
    if(mm_is_pattern(matcode))
        return MM_UNSUPPORTED_TYPE;

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
        return ret_code;

    if(mm_is_symmetric(matcode))
    {
        if (mm_is_real(matcode))
        {
            idxType i;
            for (i=0; i<nz; i++)
            {
                idxType tmpI, tmpJ;
                double tmpVal;
                if(sizeof(tmpI) == sizeof(int))
                {
                    if (fscanf(f, "%d %d %lg\n", &tmpI, &tmpJ, &tmpVal)
                        != 3) return MM_PREMATURE_EOF;
                }
                else
                if (fscanf(f, "%ld %ld %lg\n", &tmpI, &tmpJ, &tmpVal)
                    != 3) return MM_PREMATURE_EOF;

                if(tmpI != tmpJ)
                    add ++;

            }
        }
    }
    else
    {
        u_errexit("should not call from outside for nonsymmat\n");
    }
    return nz + add;

}

/************************************************************************
 mm_read_mtx_crd()  fills M, N, nz, array of values, and return
 type code, e.g. 'MCRS'

 if matrix is complex, values[] is of size 2*nz,
 (nz pairs of real/imaginary values)
 ************************************************************************/

int mm_read_mtx_crd(char *fname, idxType *M, idxType *N, idxType *nz, idxType **I, idxType **J,
                    double **val, MM_typecode *matcode)
{
    int ret_code;
    idxType bsz;

    FILE *f;

    idxType *lI;
    idxType *lJ;
    double *lval;

    idxType addind;
    if (strcmp(fname, "stdin") == 0) f=stdin;
    else
    if ((f = fopen(fname, "r")) == NULL)
        return MM_COULD_NOT_READ_FILE;


    if ((ret_code = mm_read_banner(f, matcode)) != 0)
        return ret_code;

    if (!(mm_is_valid(*matcode) && mm_is_sparse(*matcode) &&
          mm_is_matrix(*matcode)))
        return MM_UNSUPPORTED_TYPE;
    if(mm_is_complex(*matcode))
        return MM_UNSUPPORTED_TYPE;

    if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) != 0)
        return ret_code;

    if(mm_is_symmetric(*matcode))
        bsz = 2 * (*nz);
    else
        bsz = *nz;
    addind = *nz;

    lI= *I = (idxType *)  malloc(bsz * sizeof(idxType));
    lJ= *J = (idxType *)  malloc(bsz * sizeof(idxType));
    lval=*val = NULL;

#ifdef PRINT_PROGRESS
    if(sizeof(idxType) == sizeof(int))
        printf("size is read as %d, bsz=%d\n", *nz, bsz);
    else
        printf("size is read as %ld, bsz=%ld\n", *nz, bsz);
#endif
    if (mm_is_complex(*matcode))
    {
        u_errexit("did not consider the complex case\n");
    }
    else if (mm_is_real(*matcode))
    {
        lval=*val = (double *) malloc(bsz * sizeof(double));
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
                                        *matcode);
        if (ret_code != 0) return ret_code;
    }

    else if (mm_is_pattern(*matcode))
    {
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
                                        *matcode);

        if (ret_code != 0) return ret_code;
    }
    else if (mm_is_integer(*matcode))
    {
        lval=*val = (double *) malloc(bsz * sizeof(double));
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
                                        *matcode);
        if (ret_code != 0) return ret_code;
    }
    else
        return MM_UNSUPPORTED_TYPE;

    if (f != stdin) fclose(f);



    if(mm_is_symmetric(*matcode))/*add the symmetric elements*/
    {
        idxType i;
        addind = *nz;

#ifdef PRINT_PROGRESS
        printf("will add symmetric elements\n");
#endif

        for(i = 0; i < *nz; i++)
        {
            if(lI[i] != lJ[i])
            {
                lI[addind] = lJ[i];
                lJ[addind] = lI[i];
                if(!mm_is_pattern(*matcode))
                    lval[addind] = lval[i];
                addind ++;
            }
#ifdef PRINT_PROGRESS
            if(i == (*nz -1))
            {
                if(sizeof(idxType) == sizeof(int))
                    printf("addind=%d, *nz=%d, bsz=%d\n", addind, *nz, bsz);
                else
                    printf("addind=%ld, *nz=%ld, bsz=%ld\n", addind, *nz, bsz);

            }
#endif
        }
    }
    *nz = addind;

    return 0;
}

int mm_write_banner(FILE *f, MM_typecode matcode)
{
    char *str = mm_typecode_to_str(matcode);
    int ret_code;

    ret_code = fprintf(f, "%s %s\n", MatrixMarketBanner, str);
    free(str);
    if (ret_code !=2 )
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}

int mm_write_mtx_crd(char fname[], idxType M, idxType N, idxType nz, idxType I[], idxType J[],
                     double val[], MM_typecode matcode)
{
    FILE *f;
    idxType i;

    if (strcmp(fname, "stdout") == 0)
        f = stdout;
    else
    if ((f = fopen(fname, "w")) == NULL)
        return MM_COULD_NOT_WRITE_FILE;

    /* print banner followed by typecode */
    fprintf(f, "%s ", MatrixMarketBanner);
    fprintf(f, "%s\n", mm_typecode_to_str(matcode));

    /* print matrix sizes and nonzeros */
    if(sizeof(idxType) == sizeof(int))
        fprintf(f, "%d %d %d\n", M, N, nz);
    else
        fprintf(f, "%ld %ld %ld\n", M, N, nz);

    /* print values */
    if (mm_is_pattern(matcode))
    {
        for (i=0; i<nz; i++)
        {
            if(sizeof(idxType) == sizeof(int))
                fprintf(f, "%d %d\n", I[i], J[i]);
            else
                fprintf(f, "%ld %ld\n", I[i], J[i]);
        }
    }
    else if (mm_is_real(matcode))
    {
        for (i=0; i<nz; i++)
        {
            if(sizeof(idxType) == sizeof(int))
                fprintf(f, "%d %d %20.16g\n", I[i], J[i], val[i]);
            else
                fprintf(f, "%ld %ld %20.16g\n", I[i], J[i], val[i]);
        }
    }
    else if (mm_is_complex(matcode))
    {
        for (i=0; i<nz; i++)
        {
            if(sizeof(idxType) == sizeof(int))
                fprintf(f, "%d %d %20.16g %20.16g\n", I[i], J[i], val[2*i],
                        val[2*i+1]);
            else
                fprintf(f, "%ld %ld %20.16g %20.16g\n", I[i], J[i], val[2*i],
                        val[2*i+1]);
        }
    }
    else
    {
        if (f != stdout) fclose(f);
        return MM_UNSUPPORTED_TYPE;
    }

    if (f !=stdout) fclose(f);

    return 0;
}


char  *mm_typecode_to_str(MM_typecode matcode)
{
    char buffer[MM_MAX_LINE_LENGTH];
    char *types[4];
    /*    int error =0;*/
    /*    int i;*/

    /* check for MTX type */
    if (mm_is_matrix(matcode))
        types[0] = MM_MTX_STR;
    else
        return NULL;
    /* check for CRD or ARR matrix */
    if (mm_is_sparse(matcode))
        types[1] = MM_SPARSE_STR;
    else
    if (mm_is_dense(matcode))
        types[1] = MM_DENSE_STR;
    else
        return NULL;

    /* check for element data type */
    if (mm_is_real(matcode))
        types[2] = MM_REAL_STR;
    else
    if (mm_is_complex(matcode))
        types[2] = MM_COMPLEX_STR;
    else
    if (mm_is_pattern(matcode))
        types[2] = MM_PATTERN_STR;
    else
    if (mm_is_integer(matcode))
        types[2] = MM_INT_STR;
    else
        return NULL;


    /* check for symmetry type */
    if (mm_is_general(matcode))
        types[3] = MM_GENERAL_STR;
    else
    if (mm_is_symmetric(matcode))
        types[3] = MM_SYMM_STR;
    else
    if (mm_is_hermitian(matcode))
        types[3] = MM_HERM_STR;
    else
    if (mm_is_skew(matcode))
        types[3] = MM_SKEW_STR;
    else
        return NULL;

    sprintf(buffer,"%s %s %s %s", types[0], types[1], types[2], types[3]);
    return strdup(buffer);
}

int qsortfnct(const void *i1, const void *i2)
{
    idxType ii1=*(idxType*)i1, ii2 = *(idxType*)i2;
    if( ii1 < ii2) return -1;
    else if(ii1 > ii2) return 1;
    else return 0;
}

void fillCSC(idxType *xpins, idxType *pins, idxType _n, idxType _c, idxType _p, idxType *nids, idxType *cids)
{
    idxType i;
    memset(xpins, 0, sizeof(idxType)*(_n+2));

    for(i=0; i < _p; i++)
    {
        idxType nid = nids[i];
        idxType cid = cids[i];
        if(nid > _n || nid <= 0)
        {
            if(sizeof(idxType) == sizeof(int))
                u_errexit("nid=%d out of range [1...%d] from file\n", nid, _n);
            else
                u_errexit("nid=%ld out of range [1...%ld] from file\n", nid, _n);

        }
        if(cid > _c || cid <= 0)
        {
            if(sizeof(idxType) == sizeof(int))
                u_errexit("cid=%d out of range [1...%d] from file\n", cid, _c);
            else
                u_errexit("cid=%ld out of range [1...%ld] from file\n", cid, _c);

        }
        xpins[nid]++;
    }
    for(i=1; i <= _n+1; i++)
        xpins[i] += xpins[i-1];

    for(i=_p-1; i>=0; i--)
    {
        idxType nid = nids[i];
        idxType cid = cids[i];
        xpins[nid] = xpins[nid]-1;
        pins[xpins[nid]] = cid;

    }
    /*we like the indices sorted*/
    for(i=1; i<=_n; i++)
    {
        qsort(&pins[xpins[i]], xpins[i+1]-xpins[i], sizeof(idxType), qsortfnct);
    }
    if(xpins[0] != 0){
        if(sizeof(idxType) == sizeof(int))
            printf("xpins[0]=%d ~= 0\n", xpins[0]);
        else
            printf("xpins[0]=%ld ~= 0\n", xpins[0]);

        fflush(stdout);
        exit(14);
    }
    if(xpins[_n+1] != _p){
        if(sizeof(idxType) == sizeof(int))
        {
            printf("xpins[_n]=%d ~= _p=%d\n", xpins[_n+1], _p);fflush(stdout);
        }
        else
        {
            printf("xpins[_n]=%ld ~= _p=%ld\n", xpins[_n+1], _p);fflush(stdout);
        }
        exit(15);
    }
    return ;

}

void fillCSCvals(idxType *xpins, idxType *pins, ecType *vals, idxType _n, idxType _c, idxType _p, idxType *nids, idxType *cids, double *pvals)
{
    idxType i;
    memset(xpins, 0, sizeof(idxType)*(_n+2));

    for(i=0; i < _p; i++)
    {
        idxType nid = nids[i];
        idxType cid = cids[i];
        if(nid > _n || nid <= 0)
            u_errexit("nid=%d out of range [1...%d] from file\n", nid, _n);

        if(cid > _c || cid < 0)
            u_errexit("cid=%d out of range [1...%d] from file\n", cid, _c);

        xpins[nid]++;
    }
    for(i=1; i <= _n+1; i++)
        xpins[i] += xpins[i-1];

    for(i=_p-1; i>=0; i--)
    {
        idxType nid = nids[i];
        idxType cid = cids[i];
        xpins[nid] = xpins[nid]-1;
        pins[xpins[nid]] = cid;
        if (pvals)
            vals[xpins[nid]] = (ecType)pvals[i];
    }

    if(xpins[1] != 0){
        printf("xpins[0]=%d ~= 0\n", xpins[1]);
        fflush(stdout);
        exit(14);
    }
    if(xpins[_n+1] != _p){
        printf("xpins[_n]=%d ~= _p=%d\n", xpins[_n+1], _p);fflush(stdout);
        exit(15);
    }
    return ;

}

int mm_read_into_csc_pattern(char *fname, idxType *M, idxType *N, idxType *nnz, idxType **col_ptrs, idxType **row_inds)
{
    idxType *rindx=NULL, *cindx=NULL;
    double *vals=NULL;
    int  rcode;
    MM_typecode matcode;
    rcode = mm_read_mtx_crd(fname, M, N, nnz, &rindx, &cindx,
                            &vals, &matcode);
    if(rcode != 0)
        return rcode;
#ifdef PRINT_PROGRESS
    if(sizeof(idxType) == sizeof(int))
        printf("bmmio: will fill in csc storage arrays %d %d %d\n", *M, *N, *nnz);
    else
        printf("bmmio: will fill in csc storage arrays %ld %ld %ld\n", *M, *N, *nnz);

#endif

    *col_ptrs= (idxType *) malloc(sizeof(idxType)*(*N+2));
    *row_inds = (idxType *) malloc(sizeof(idxType)*(*nnz));
    fillCSC(*col_ptrs, *row_inds, *N, *M, *nnz, cindx, rindx);
    free(rindx);
    free(cindx);
    free(vals);
    return 0;
}

int mm_read_into_csc(char *fname, idxType *M, idxType *N, idxType *nnz, idxType **col_ptrs, idxType **row_inds, ecType **vals)
{
    idxType *rindx=NULL, *cindx=NULL;
    double *lvals=NULL;
    int  rcode;
    MM_typecode matcode;
    rcode = mm_read_mtx_crd(fname, M, N, nnz, &rindx, &cindx,
                            &lvals, &matcode);
    if(rcode != 0)
        return rcode;


#ifdef PRINT_PROGRESS
    if(sizeof(idxType) == sizeof(int))
        printf("bmmio: will fill in csc storage arrays %d %d %d\n", *M, *N, *nnz);
    else
        printf("bmmio: will fill in csc storage arrays %ld %ld %ld\n", *M, *N, *nnz);

#endif

    *col_ptrs= (idxType *) malloc(sizeof(idxType)*(*N+2));
    *row_inds = (idxType *) malloc(sizeof(idxType)*(*nnz));
    if (mm_is_pattern(matcode))
        *vals = NULL;
    else
        *vals = (ecType *) malloc(sizeof(ecType) * (*nnz));


    fillCSCvals(*col_ptrs, *row_inds, *vals, *N, *M, *nnz, cindx, rindx, lvals);
    free(rindx);
    free(cindx);
    free(lvals);
    return 0;
}
#pragma GCC diagnostic pop


int loadFromCSC(dgraph *G, idxType nVrtx, idxType nEdge, idxType *col_ptrs, idxType *row_inds, ecType *vals)
{
    /* col_ptrs [1...nVrtx+1] are used.
     * row_inds [0..nEdges-1] are used and contains indices 1...n (numbering are fortran like).
     *
     * assumes G's data allocated and G->nEdge > nEdge.
     *
     * at the end we set G->nEdge to the exact number of edges, where the lower traingular,
     * including the diagonal, entries are discarded from the CSC structure.
     *
     * the edgeCosts (in vals) are mode positive, if any.
     *
     * (i,j) is an edge from i to j. Therefore the col_ptrs[i] and row_inds[colptrs[i]] to
     * colptrs[i+1]-1] define incoming edges of i.
     */
    idxType i, at, j;
    printf ("\n---->%d %d %d %d\n", G->nVrtx, nVrtx, G->nEdge, nEdge);

    if (G->nVrtx < nVrtx || G->nEdge < nEdge) {
        printf ("\n---->%d %d %d %d\n", G->nVrtx, nVrtx, G->nEdge, nEdge);
        u_errexit("loadFromCSC: potentially not enough space in G");
    }

    G->nVrtx = nVrtx;
    G->inStart[0] = 0;
    G->inStart[1] = at = 0;
    long upctr=0,lowctr=0;
    /*inStart[i]  is already set*/
    for (j = 1; j<=nVrtx; j++)
    {
        for (i = col_ptrs[j]; i < col_ptrs[j+1]; i++)
        {
            if (row_inds[i] < j)
            {
                ++upctr;
            }
            if (row_inds[i]>j) {
                ++lowctr;
            }
        }
    }
    if (upctr>=lowctr) {
        /*inStart[i]  is already set*/
        for (j = 1; j<=nVrtx; j++) {
            for (i = col_ptrs[j]; i < col_ptrs[j+1]; i++) {
                if (row_inds[i] < j) {
                    G->in[at++] = row_inds[i];
                    if ( (G->frmt == 2 || G->frmt == 3) && (vals))
                        G->ecIn[at-1] = vals[i] >=0 ? vals[i]: -vals[i];
                    else
                        G->ecIn[at-1] = 1;
                }
            }
            G->inEnd[j] = at-1;
            G->inStart[j+1] = at;
        }
    }
    else{
        /*inStart[i]  is already set*/
        for (j = 1; j<=nVrtx; j++) {
            for (i = col_ptrs[j]; i < col_ptrs[j+1]; i++) {
                if (row_inds[i] > j) {
                    G->in[at++] = row_inds[i];
                    if ( (G->frmt == 2 || G->frmt == 3) && (vals))
                        G->ecIn[at-1] = vals[i] >=0 ? vals[i]: -vals[i];
                    else
                        G->ecIn[at-1] = 1;
                }
            }
            G->inEnd[j] = at-1;
            G->inStart[j+1] = at;
        }
    }
    G->nEdge = at;

    G->nbsources = 0;
    G->nbtargets = 0;

    fillOutFromIn(G);
    for (i=1; i<=nVrtx; i++){
        if (G->inEnd[i] < G->inStart[i])
            G->sources[G->nbsources++] = i;
        if (G->outEnd[i] < G->outStart[i])
            G->targets[G->nbtargets++] = i;
    }
    sortNeighborLists(G);
    return 0;
}

int readDGraph(dgraph *G, char* file_name, int prefer_binary)
/* readDGraph will call following four based on the extensions:
    ".dg" or ".dgraph", ".mtx", ".bin", ".dot"
    and if prefer_binary is one, it will read from binary, and generate binary if it wasn't binary */
{
    int  ret=0;
    char fnbin[PATH_MAX];

    strcpy(fnbin, file_name);
    strcat(fnbin, ".bin");

    if (strstr(file_name, ".bin")) { /* file is already binary */
        return readDGraphBinary(G, file_name);
    } else {

        if (prefer_binary) {
            FILE *fp = fopen(fnbin, "rb");

            if (fp!=NULL) {
                fclose(fp);
                return readDGraphBinary(G, fnbin);
            }
        }
        if (strstr(file_name, ".dot")) {
            ret = readDGraphDot(G, file_name);

        }else if (strstr(file_name, ".mtx")) {
            ret = readDGraphMtx(G, file_name);
        }else if (strstr(file_name, ".fab")) {
            ret = readDGraphFab(G, file_name);
        }else /* no extension or .dg or .dgraph */{
            ret = readDGraphDG(G, file_name);
        }
    }

    if (!ret && prefer_binary)
        ret = writeDGraphBinary(G, fnbin);
    return ret;
}

int readDGraphDG(dgraph *G, char* file_name)
{
    FILE* file;
    file = fopen(file_name, "r");
    if (file == NULL)
        u_errexit("readDGraphDG: cannot open %s\n",file_name);
    int nVrtx, nEdge, i, tmp;
    fscanf(file, "%d", &nVrtx);
    fscanf(file, "%d", &nEdge);

    G->frmt = 3;
    G->nVrtx = nVrtx;
    G->nEdge = nEdge;
    G->inStart  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->hollow  = (int * ) calloc((nVrtx + 2),sizeof(int));
    G->inEnd= (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->in = (idxType * ) malloc(sizeof(idxType) * nEdge);
    G->totvw = nVrtx;
    G->maxindegree = 0;

    G->inStart[0] = 0;
    for (i=1; i<=nVrtx+1; i++) {
        fscanf(file, "%d", &tmp);
        G->inStart[i] = tmp;
    }

    G->inEnd[0] = -1;
    for (i=1; i<=nVrtx; i++) {
        G->inEnd[i] = G->inStart[i+1] -1;
        int degree = G->inEnd[i] - G->inStart[i] + 1;
        G->maxindegree = G->maxindegree < degree ? degree : G->maxindegree;
    }
    G->inEnd[nVrtx+1] = nEdge;
    for (i=0; i<nEdge; i++) {
        fscanf(file, "%d", &tmp);
        G->in[i] = tmp;
    }

    G->vw = NULL;
    G->ecIn = G->ecOut = NULL;

    if (nEdge <= 0)
        u_errexit("allocateDGraphData: empty edge set\n");

    if (nVrtx <=-1)
        u_errexit("allocateDGraphData: empty vertex set\n");

    G->outStart  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->outEnd  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->out = (idxType * ) malloc(sizeof(idxType) * nEdge);

    G->vw = (vwType *) malloc(sizeof(vwType) * (nVrtx+1));
    for (i=1; i<=nVrtx; i++)
        G->vw[i] = 1;
    G->ecIn = (ecType *) malloc(sizeof(ecType) * nEdge);
    G->ecOut = (ecType *) malloc(sizeof(ecType) * nEdge);
    for (i=0; i< nEdge; i++) {
        G->ecIn[i] = 1;
        G->ecOut[i] = 1;
    }

    fillOutFromIn(G);

    G->sources  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 1));
    G->targets  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 1));
    G->nbsources = 0;
    G->nbtargets = 0;
    for (i=1; i<=nVrtx; i++){
        if (G->inEnd[i] < G->inStart[i])
            G->sources[G->nbsources++] = i;
        if (G->outEnd[i] < G->outStart[i])
            G->targets[G->nbtargets++] = i;
    }

    G->maxVW = G->maxEC = -1;
    fclose(file);
    return 0;
}

int readDGraphMtx(dgraph *G, char* file_name)
{
    idxType i;
    int rcodeMM;

    idxType M, N, nnz, *col_ptrs, *row_inds;
    M=N=nnz=-1;
    rcodeMM =  mm_read_into_csc_pattern(file_name, &M, &N, &nnz, &col_ptrs, &row_inds);
    if (rcodeMM != 0) {
        u_errexit("Cannot read file %s\n", file_name);
        return -1;
    }

    if (M != N)
        u_errexit("%s is not a square matrix\n", file_name);
    allocateDGraphData(G,M,nnz,3);
    /*the vertex weights are not set, while 3 is used*/
    /*let us set them*/

    if (G->vw == NULL) {
        G->vw = (vwType*) malloc(sizeof(vwType)*(M+1));
    }
    for (i=0; i <= M; i++)
        G->vw[i] = 1;

    return loadFromCSC(G, M, nnz, col_ptrs, row_inds, NULL);
}


int readDGraphBinary(dgraph *G, char* file_name)
{
    FILE *fp = fopen(file_name,"rb");

    if (fp == NULL) {
        u_errexit("Cannot open file %s\n", file_name);
    }

    fread(&G->frmt, sizeof(G->frmt), 1, fp);
    fread(&G->nVrtx, sizeof(G->nVrtx), 1, fp);
    fread(&G->nEdge, sizeof(G->nEdge), 1, fp);
    fread(&G->totvw, sizeof(G->totvw), 1, fp);
    fread(&G->maxindegree, sizeof(G->maxindegree), 1, fp);
    fread(&G->maxoutdegree, sizeof(G->maxoutdegree), 1, fp);

    idxType nVrtx = G->nVrtx;
    idxType nEdge = G->nEdge;

    if (nEdge <= 0)
        u_errexit("allocateDGraphData: empty edge set\n");

    if (nVrtx <=-1)
        u_errexit("allocateDGraphData: empty vertex set\n");

    G->hollow  = (int * ) calloc(nVrtx + 2, sizeof(int));
    G->inStart  = (idxType * ) calloc(nVrtx + 2, sizeof(idxType));
    G->inEnd= (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->in = (idxType * ) malloc(sizeof(idxType) * nEdge);

    G->outStart  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->outEnd  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->out = (idxType * ) malloc(sizeof(idxType) * nEdge);

    G->sources  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 1));
    G->targets  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 1));

    G->vw = (vwType *) malloc(sizeof(vwType) * (nVrtx+1));

    G->ecIn = (ecType *) malloc(sizeof(ecType) * nEdge);
    G->ecOut = (ecType *) malloc(sizeof(ecType) * nEdge);

    fread(G->inStart, sizeof(G->inStart[0]), (nVrtx+2), fp);
    fread(G->inEnd, sizeof(G->inEnd[0]), (nVrtx+2), fp);
    fread(G->in, sizeof(G->in[0]), nEdge, fp);
    fread(G->vw, sizeof(G->vw[0]), (nVrtx+1), fp);
    fread(G->outStart, sizeof(G->outStart[0]), (nVrtx+2), fp);
    fread(G->outEnd, sizeof(G->outEnd[0]), (nVrtx+2), fp);
    fread(G->out, sizeof(G->out[0]), nEdge, fp);
    fread(G->ecIn, sizeof(ecType), nEdge, fp);
    fread(G->ecOut, sizeof(ecType), nEdge, fp);

    fread(&G->nbsources, sizeof(G->nbsources), 1, fp);
    fread(&G->nbtargets, sizeof(G->nbtargets), 1, fp);
    fread(G->sources, sizeof(G->sources[0]), nVrtx, fp);
    fread(G->targets, sizeof(G->targets[0]), nVrtx, fp);

    fclose(fp);
    fillOutFromIn(G);

    return 0;
}

int writeDGraphBinary(dgraph *G, char* file_name)
{
    FILE *fp = fopen(file_name, "wb");

    if (fp == NULL) {
        u_errexit("writeDGraphBinary: Cannot open file %s for writing\n", file_name);
        return -1;
    }

    fwrite(&G->frmt, sizeof(G->frmt), 1, fp);
    fwrite(&G->nVrtx, sizeof(G->nVrtx), 1, fp);
    fwrite(&G->nEdge, sizeof(G->nEdge), 1, fp);
    fwrite(&G->totvw, sizeof(G->totvw), 1, fp);
    fwrite(&G->maxindegree, sizeof(G->maxindegree), 1, fp);
    fwrite(&G->maxoutdegree, sizeof(G->maxoutdegree), 1, fp);
    fwrite(G->inStart, sizeof(G->inStart[0]), (G->nVrtx+2), fp);
    fwrite(G->inEnd, sizeof(G->inEnd[0]), (G->nVrtx+2), fp);
    fwrite(G->in, sizeof(G->in[0]), G->nEdge, fp);
    fwrite(G->vw, sizeof(G->vw[0]), (G->nVrtx+1), fp);
    fwrite(G->outStart, sizeof(G->outStart[0]), (G->nVrtx+2), fp);
    fwrite(G->outEnd, sizeof(G->outEnd[0]), (G->nVrtx+2), fp);
    fwrite(G->out, sizeof(G->out[0]), G->nEdge, fp);
    fwrite(G->ecIn, sizeof(ecType), G->nEdge, fp);
    fwrite(G->ecOut, sizeof(ecType), G->nEdge, fp);
    fwrite(&G->nbsources, sizeof(G->nbsources), 1, fp);
    fwrite(&G->nbtargets, sizeof(G->nbtargets), 1, fp);
    fwrite(G->sources, sizeof(G->sources[0]), (G->nVrtx+1), fp);
    fwrite(G->targets, sizeof(G->targets[0]), (G->nVrtx+1), fp);

    fclose(fp);
    return 0;
}

int writePartsToTxt(dgraph *G, char* file_name, idxType* part){
    FILE *file;
    file = fopen(file_name, "w");
    idxType i;
    for (i=1; i<=G->nVrtx; i++) {
        fprintf(file, "%d\n", (int) part[i]);
    }
    fclose(file);
    return 0;
}

/* part can be NULL, if not null colors generated for each part */
int writeDGraphDot(dgraph *G, char* file_name, idxType* part)
{
    FILE    *file;
    idxType i, j;
    char*   col[30];
    col[0] = "red";
    col[1] = "blue";
    col[2] = "green";
    col[3] = "purple";
    col[4] = "orange";
    col[5] = "brown";
    col[6] = "chocolate";
    col[7] = "crimson";
    col[8] = "salmon";
    col[9] = "grey";
    col[10] = "gold";
    col[11] = "indigo";

    file = fopen(file_name, "w");
    if (file == NULL)
        u_errexit("writeDGraphDot: cannot open %s\n", file_name);
    fprintf(file, "digraph G {\n");
    for (i=1; i<=G->nVrtx; i++) {
        if (part)
            fprintf(file, "%d [style=filled, fillcolor=%s, weight=%d];\n", i-1, col[part[i] % 12], G->vw[i]);
        else
            fprintf(file, "%d [weight=%d];\n", i-1, G->vw[i]);
    }
    for (i=1; i<=G->nVrtx; i++) {
        for (j=G->outStart[i]; j<=G->outEnd[i]; j++) {
            idxType outnode = G->out[j];
            fprintf(file, "%d->%d [weight=%d];\n", i-1, outnode-1, (int)G->ecOut[j]);
        }
    }
    fprintf(file, "}\n");
    fclose(file);

    return 0;
}

void writePartGraphDot(dgraph* G, char* file_name, idxType* part)
{
    FILE    *file;
    int maxParts= -1, i = 1, j=0;
    for (i=1; i<=G->nVrtx; ++i)
        if (maxParts< part[i])
            maxParts = part[i];
    int** adj = (int**) calloc(maxParts+1, sizeof(int*));
    for (i=0; i<maxParts; ++i)
        adj[i] = (int*) calloc(maxParts+1,sizeof(int));
    for (i=1; i<=G->nVrtx; ++i)
        for(j=G->outStart[i]; j<=G->outEnd[i]; ++j) {
            if (part[i] != part[G->out[j]])
                ++ adj[part[i]][part[G->out[j]]];
            // if (part[i] == 77 || part[i]==78)
            //     printf("node %d part %d out %d part %d\n", i, part[i], G->out[j], part[G->out[j]]);
        }

    file = fopen(file_name, "w");
    if (file == NULL)
        u_errexit("writePartGraphDot: cannot open %s\n", file_name);
    fprintf(file, "digraph G {\n");
    for (i=0; i<maxParts; i++) {
        fprintf(file, "%d;\n", i);
    }
    for (i=0; i<maxParts; i++) {
        for (j=0; j<maxParts; j++) {
            if(adj[i][j] > 0)
                fprintf(file, "%d->%d [weight=%d];\n", i, j, adj[i][j]);
        }
    }
    fprintf(file, "}\n");
    fclose(file);
    for (i=0; i<maxParts; ++i)
        free(adj[i]);
    free(adj);
    return;

}

/* part can be NULL, if not null colors generated for each part */
int writeDGraphDotWithProc(dgraph *G, char* file_name, idxType* part, idxType* proc)
{
    FILE    *file;
    idxType i, j;
    char*   col[30];
    col[0] = "red";
    col[1] = "blue";
    col[2] = "green";
    col[3] = "purple";
    col[4] = "orange";
    col[5] = "brown";
    col[6] = "chocolate";
    col[7] = "crimson";
    col[8] = "salmon";
    col[9] = "grey";
    col[10] = "gold";
    col[11] = "indigo";

    file = fopen(file_name, "w");
    if (file == NULL)
        u_errexit("writeDGraphDot: cannot open %s\n", writeDGraphDot);
    fprintf(file, "digraph G {\n");
    for (i=1; i<=G->nVrtx; i++) {
        if (part!=NULL && proc!=NULL)
            fprintf(file, "%d [style=filled, penwidth=8, color=%s, fillcolor=%s, weight=%d];\n", i-1, col[part[i] % 12], col[proc[i] % 12],G->vw[i]);
        else
            fprintf(file, "%d [weight=%d];\n", i-1, G->vw[i]);
    }
    for (i=1; i<=G->nVrtx; i++) {
        for (j=G->outStart[i]; j<=G->outEnd[i]; j++) {
            idxType outnode = G->out[j];
            fprintf(file, "%d->%d [weight=%d];\n", i-1, outnode-1, (int)G->ecOut[j]);
        }
    }
    fprintf(file, "}\n");
    fclose(file);

    return 0;
}


#define DOT_LINE_MAX 1000

/* note that this ignores whole like if it has any of these chars */
static inline int dotCommentLine(char *line)
{
    return (line[0] == '#') || strstr(line, "//") || strstr(line, "/*") || strstr(line, "*/");
}



int valid_integer(char c){
    if ((c != '0') && (c != '1') && (c != '2') && (c != '3') && (c != '4') && (c != '5') && (c != '6') && (c != '7') && (c != '8') && (c != '9'))
        return 0;
    return 1;
}

int readDGraphFab(dgraph *G, char* file_name)
{
    FILE* file = fopen(file_name, "r");
    char delims[] = " ->[,=];\n";
    int graphStatus=0; /* 0 not started, 1 saw '{' */

    if (file == NULL)
        u_errexit("readDGraphFab: Cannot open file %s\n", file_name);

    char line[DOT_LINE_MAX];
    idxType nVrtx = 0, nEdge = 0, innode, outnode;
    idxType i, j;
    idxType* nbneighbors=NULL;
    idxType** inneighbors=NULL;
    char* pos;
    char* token;

    while (!feof(file)) {
        do {
            fgets(line, DOT_LINE_MAX, file);
        } while (!feof(file) && dotCommentLine(line));
        if (feof(file))
            break;

        if (strlen(line) >= LINE_MAX - 2)
            u_errexit("readDGraphFab: line is too long! \n'%'\n", line);

        /* it is an edge */
        pos = strstr(line, ":");
        if (pos == NULL) {
            continue;
        }

        token = strtok(line, ":");
        if (valid_integer(token[0]) == 0)
            u_errexit("readDGraphFab: first node not valid! \n'%'\n", line);

        innode = atoi(token);
        nVrtx = nVrtx < innode ? innode + 1 : nVrtx;
        token = strtok(line, " ");
        while (token != NULL) {
            if (token[0] == '#')
                break;
            if (valid_integer(token[0]) == 0)
                u_errexit("readDGraphFab: output node not valid! \n'%'\n", line);
            outnode = atoi(token);
            nEdge++;
            nVrtx = nVrtx < outnode ? outnode + 1 : nVrtx;
            token = strtok(line, " ");
        }
    }
    fclose(file);

    nbneighbors = (idxType*) calloc(nVrtx+1, sizeof(idxType));
    inneighbors = (idxType**) malloc((nVrtx+1) * sizeof(idxType*));
    for (i = 1; i <= nVrtx; i++)
        inneighbors[i] = (idxType*) malloc(100 * sizeof(idxType));
    file = fopen(file_name, "r");

    while (!feof(file)) {
        do {
            fgets(line, DOT_LINE_MAX, file);
        } while (!feof(file) && dotCommentLine(line));
        if (feof(file))
            break;

        pos = strstr(line, ":");
        if (pos == NULL) {
            continue;
        }

        char *token = strtok(line, ":");
        innode = atoi(token);
        token = strtok(line, " ");
        while (token != NULL) {
            outnode = atoi(token);
            token = strtok(line, " ");
        }
        if ((nbneighbors[outnode + 1] + 1) % 100 == 0)
            inneighbors[outnode + 1] = (idxType *) realloc(inneighbors[outnode + 1],
                                                           (nbneighbors[outnode + 1] + 101) * sizeof(idxType));
        inneighbors[outnode + 1][nbneighbors[outnode + 1]++] = innode + 1;
    }
    fclose(file);

    G->frmt = 3;
    G->nVrtx = nVrtx;
    G->nEdge = nEdge;
    G->totvw = nVrtx;
    G->maxindegree = 0;
    G->hollow  = (int * ) calloc(nVrtx + 2, sizeof(int));
    G->inStart  = (idxType * ) calloc(nVrtx + 2, sizeof(idxType));
    G->inEnd= (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->in = (idxType * ) malloc(sizeof(idxType) * nEdge);
    idxType idx = 0, degree;

    for (i=1; i<=nVrtx; i++) {
        G->inStart[i] = idx;
        G->inEnd[i-1] = idx-1;
        if (i>1) {
            degree = G->inEnd[i-1] - G->inStart[i-1] + 1;
            G->maxindegree = G->maxindegree < degree ? degree : G->maxindegree;
        }
        for (j=0; j< nbneighbors[i]; j++) {
            G->in[idx++] = inneighbors[i][j];
        }
    }
    G->inStart[0] = 0;
    G->inEnd[0] = -1;
    G->inEnd[nVrtx] = idx-1;
    G->inEnd[nVrtx+1] = nEdge;
    G->inStart[nVrtx+1] = nEdge;
    degree = G->inEnd[nVrtx] - G->inStart[nVrtx] + 1;
    G->maxindegree = G->maxindegree < degree ? degree : G->maxindegree;
    if (nEdge <= 0)
        u_errexit("allocateDGraphData: empty edge set\n");

    if (nVrtx <=-1)
        u_errexit("allocateDGraphData: empty vertex set\n");

    G->outStart  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->outEnd  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->out = (idxType * ) malloc(sizeof(idxType) * nEdge);

    G->vw = (vwType *) malloc(sizeof(vwType) * (nVrtx+1));
    for (i=1; i<=nVrtx; i++)
        G->vw[i] = 1;
    G->ecIn = (ecType *) malloc(sizeof(ecType) * nEdge);
    G->ecOut = (ecType *) malloc(sizeof(ecType) * nEdge);
    for (i=0; i< nEdge; i++) {
        G->ecIn[i] = 1;
        G->ecOut[i] = 1;
    }

    fillOutFromIn(G);

    G->sources  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 1));
    G->targets  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 1));
    G->nbsources = 0;
    G->nbtargets = 0;
    for (i=1; i<=nVrtx; i++){
        if (G->inEnd[i] < G->inStart[i])
            G->sources[G->nbsources++] = i;
        if (G->outEnd[i] < G->outStart[i])
            G->targets[G->nbtargets++] = i;
    }

    for (i = 1; i <= nVrtx; i++)
        free(inneighbors[i]);
    free(nbneighbors);
    free(inneighbors);
    return 0;
}

idxType readVector(char* file_name, idxType* vector)
{
    FILE* file = fopen(file_name, "r");
    char delims[] = " ->[,=];\n";

    if (file == NULL)
        u_errexit("readVector: Cannot open file %s\n", file_name);

    char line[DOT_LINE_MAX];
    idxType node;
    idxType i, j;
    char* token;
    idxType nbvector = 1;

    while (!feof(file)) {
        do {
            fgets(line, DOT_LINE_MAX, file);
        } while (!feof(file) && dotCommentLine(line));
        if (feof(file))
            break;

        if (strlen(line) >= LINE_MAX - 2)
            u_errexit("readVector: line is too long! \n'%'\n", line);

        /* it is an edge */
        node = atoi(line);
        vector[nbvector++] = node;
    }
    fclose(file);

    return nbvector-1;
}
