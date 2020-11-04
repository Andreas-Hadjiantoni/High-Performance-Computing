#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"

// Define output file name
#define OUTPUT_FILE "stencil.pgm"

void stencilSyncGrid(const int nx, const int ny, float *  image, float *  tmp_image, int width, int height, float* halosSend, float* halosRecv, int edgesTouching[4]);
void stencilSyncCols(const int nx, const int ny, float *  image, float *  tmp_image, int width);
void stencilAsyncCols(const int nx, const int ny, float * restrict image, float * restrict tmp_image, int width, MPI_Request* requests);
void init_image(const int nx, const int ny, float *  image, float *  tmp_image);
void output_image(const char * file_name, const int nx, const int ny, float *image);
double wtime(void);

int rank;               /* 'rank' of process among it's cohort */
int size;               /* size of cohort, i.e. num processes started */
int flag;               /* for checking whether MPI_Init() has been called */
enum bool  {FALSE, TRUE};
int strlength;             /* length of a character array */
char hostname[MPI_MAX_PROCESSOR_NAME]; /* character array to hold hostname running process */
MPI_Comm comm_cart; /* a cartesian topology aware communicator */

int main(int argc, char *argv[]) {

  /*

  The code can be run in three modes:
  1) Synchronous data exchnge
  2) Asynchronous data exchnge
  3) Synchronous data Exchange with cartesian Topology

  To choose a mode, set the variable "runningMode" to the values 1, 2, or 3. Each
  of these values will correspond to the respective mode in the ordered outlined above.
  */

  int runningMode = 1;

  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  MPI_Init(&argc, &argv);

  MPI_Initialized(&flag);
  if ( flag != TRUE ) {
  	printf("failed");
  	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }

  MPI_Get_processor_name(hostname,&strlength);

  MPI_Comm_size( MPI_COMM_WORLD, &size );

  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Initiliase problem dimensions from command line arguments
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int niters = atoi(argv[3]);

  float *image;
  float *tmp_image;

  if (rank == 0) {
    // Allocate the image
    image = malloc(sizeof(float)*nx*ny);
    tmp_image = malloc(sizeof(float)*nx*ny);

    // Set the input image
    init_image(nx, ny, image, tmp_image);
  }

  int grid;
  int width;
  int height;

  float* chunk;
  float* tmp_chunk;

  float *halosSend;
  float *halosRecv;

  MPI_Request requests[4];

  int edgesTouching[4] = {0, 0, 0, 0};
  int dims[2];       /* array to hold dimensions of an NDIMS grid of processes */


  if (runningMode == 3)
    grid = TRUE;
  else
    grid = FALSE;

  if (grid) {

      int periods[2];    /* array to specificy periodic boundary conditions on each dimension */
      int coords[2];     /* array to hold the grid coordinates for a rank */

      for (int i = 0; i < 2; i++) {
        periods[i] = 0;
        dims[i] = 0;
      }

      MPI_Dims_create(size, 2, dims);

      MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_cart);

      MPI_Cart_coords(comm_cart, rank, 2, coords);

      width = nx / dims[1];
      height = ny / dims[0];

      /*
      edgesTouching[0]: 1 if rank touches top edge, 0 otherwise
      edgesTouching[1]: 1 if rank touches bottom edge, 0 otherwise
      edgesTouching[2]: 1 if rank touches left edge, 0 otherwise
      edgesTouching[3]: 1 if rank touches right edge, 0 otherwise
      */

      if (coords[0] == 0) edgesTouching[0] = 1;
      if (coords[1] == 0) edgesTouching[2] = 1;
      if (coords[0] == dims[0] - 1) edgesTouching[1] = 1;
      if (coords[1] == dims[1] - 1) edgesTouching[3] = 1;

      chunk = malloc(sizeof(float) * width * height);
      tmp_chunk = malloc(sizeof(float) * width * height);

      halosSend = malloc(sizeof(float) * 2 * (width + height));
      halosRecv = malloc(sizeof(float) * 2 * (width + height));


      if (rank == 0) {
        for (int i = 0; i < nx; i++) {
          for (int h = 0; h < dims[0]; h++) {
            //skip sending to rank 0
            if (i < width && h == 0) {
              memcpy(chunk + height * i, image + i * ny, sizeof(float) * height);
              h = 1;
            }

            int rankToSend;
            int sendCoords[2] = {h, i / width};

            MPI_Cart_rank(comm_cart, sendCoords, &rankToSend);
            MPI_Send(image + i * ny + h * height, height, MPI_FLOAT, rankToSend, 0, comm_cart);
          }

        }

      } else {

        for (int i = 0; i < width; i++) {
          MPI_Recv(chunk + height * i, height, MPI_FLOAT, 0, 0, comm_cart, MPI_STATUS_IGNORE);
        }

      }

  } else {

    // Columns:
    width = nx / size;

    int sendCounts[size];
    int displs[size];

    int widthToBeIncreased = FALSE;

    int rem = nx % size;

    if (rem > rank)
      widthToBeIncreased = TRUE;

    for (int i = 0; i < size; i++) {
      if (rem > 0) {
        rem--;
        sendCounts[i] = (width + 1) * ny;
      } else
      sendCounts[i] = width * ny;

      if (i == 0)
      displs[i] = 0;
      else
      displs[i] = displs[i - 1] + sendCounts[i - 1];
    }

    if (widthToBeIncreased == TRUE)
      width++;

    chunk = malloc(sizeof(float) * (width + 2) * ny);
    tmp_chunk = malloc(sizeof(float) * (width + 2) * ny);

    MPI_Scatterv(image, sendCounts, displs, MPI_FLOAT, chunk + ny, width * ny, MPI_FLOAT, 0, MPI_COMM_WORLD);

  }

  double tic = wtime();

  switch (runningMode) {
    case 1: {
      for (int t = 0; t < niters; ++t) {
        stencilSyncCols(nx, ny, chunk, tmp_chunk, width);
      	stencilSyncCols(nx, ny, tmp_chunk, chunk, width);
    	}
      break;
    }
    case 2: {
      for (int t = 0; t < niters; ++t) {
        stencilAsyncCols(nx, ny, chunk, tmp_chunk, width, requests);
      	stencilAsyncCols(nx, ny, tmp_chunk, chunk, width, requests);
    	}
      break;
    }
    case 3: {
      for (int t = 0; t < niters; ++t) {
        stencilSyncGrid(nx, ny, chunk, tmp_chunk, width, height, halosSend, halosRecv, edgesTouching);
        stencilSyncGrid(nx, ny, tmp_chunk, chunk, width, height, halosSend, halosRecv, edgesTouching);
      }
      break;
    }
    default:
      break;
  }


  double toc = wtime();

  if (grid == TRUE) {
    if (rank == 0) {
      for (int i = 0; i < nx; i++)
        for (int h = 0; h < dims[0]; h++) {
          //skip receiving from rank 0
          if (i < width && h == 0) {
            memcpy(image + i * ny, chunk + height * i, sizeof(float) * height);
            h = 1;
          }

          int rankToRecv;
          int recvCoords[2] = {h, i / width};

          MPI_Cart_rank(comm_cart, recvCoords, &rankToRecv);
          MPI_Recv(image + i * ny + h * height, height, MPI_FLOAT, rankToRecv, 0, comm_cart, MPI_STATUS_IGNORE);
        }
    } else
    for (int i = 0; i < width; i++)
      MPI_Send(chunk + height * i, height, MPI_FLOAT, 0, 0, comm_cart);


  } else
    MPI_Gather(chunk + ny, ny * width, MPI_FLOAT, image, ny * width, MPI_FLOAT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    // Output
    printf("------------------------------------\n");
    printf("runtime: %lf s\n", toc-tic);
    printf("------------------------------------\n");

    output_image(OUTPUT_FILE, nx, ny, image);

    free(image);
    free(tmp_image);
  }

  MPI_Finalize();

}

void stencilSyncGrid(
  const int nx,
  const int ny,
  float * restrict image,
  float * restrict tmp_image,
  int width,
  int height,
  float* halosSend,
  float* halosRecv,
  int edgesTouching[4]) {

  /*
  edgesTouching[0]: 1 if rank touches top edge, 0 otherwise
  edgesTouching[1]: 1 if rank touches bottom edge, 0 otherwise
  edgesTouching[2]: 1 if rank touches left edge, 0 otherwise
  edgesTouching[3]: 1 if rank touches right edge, 0 otherwise
  */
  for (int i = 0; i < width; i++) {
    halosSend[i] = image[i * height];
    halosSend[i + width] = image[(i + 1) * height -1];
  }

  memcpy(halosSend + 2 * width, image, sizeof(float) * height);
  memcpy(halosSend + 2 * width + height, image + (width - 1) * height, sizeof(float) * height);

  int counts[4] = {width, width, height, height};
  int displs[4] = {0, width, 2*width, 2*width + height};

  //sendrecv here:
  //MPI_Neighbor_alltoall(halosSend, width, MPI_FLOAT, halosRecv, width, MPI_FLOAT, comm_cart);
  MPI_Neighbor_alltoallv(halosSend, counts, displs, MPI_FLOAT, halosRecv, counts, displs, MPI_FLOAT, comm_cart);
  // 2 * (width + height)

  float* haloTop = halosRecv;
  float* haloBottom = halosRecv + width;
  float* haloLeft = haloBottom + width;
  float* haloRight = haloLeft + height;

  for (int i = 0; i < width; ++i) {
    for (int j = 0; j < height; ++j) {
      tmp_image[j+i*height] = image[j+i*height] * 0.6f;

      if (i > 0)
        tmp_image[j+i*height] += image[j  +(i-1)*height] * 0.1f;
      else
        if (edgesTouching[2] == 0)
          tmp_image[j + i * height] += haloLeft[j] * 0.1f;

      if (i < width - 1)
        tmp_image[j+i*height] += image[j  +(i+1)*height] * 0.1f;
      else
        if (edgesTouching[3] == 0)
          tmp_image[j+i*height] += haloRight[j] * 0.1f;

      if (j > 0)
        tmp_image[j+i*height] += image[j-1+i*height] * 0.1f;
      else
        if (edgesTouching[0] == 0)
          tmp_image[j+i*height] += haloTop[i] * 0.1f;

      if (j < height - 1)
        tmp_image[j+i*height] += image[j+1+i*height] * 0.1f;
      else
        if (edgesTouching[1] == 0)
          tmp_image[j+i*height] += haloBottom[i] * 0.1f;

    }
  }
}

void stencilAsyncCols(const int nx, const int ny, float * restrict image, float * restrict tmp_image, int width, MPI_Request *requests) {

    /*
      requests[0]: SendLeft
      requests[1]: RecvLeft
      requests[2]: sendRight
      requests[3]: RecvRight
    */
    if (rank % 2 == 1) {
      if (rank < size - 1) {
        //send left
        MPI_Isend(image + ny, ny, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &requests[0]);
        //recv left
        MPI_Irecv(image, ny, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &requests[1]);

        //send right
        MPI_Isend(image + (width * ny), ny, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &requests[2]);
        //recv right
        MPI_Irecv(image + ((width + 1) * ny), ny, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &requests[3]);
      } else {
        //send left
        MPI_Isend(image + ny, ny, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &requests[0]);
        //recv left
        MPI_Irecv(image, ny, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &requests[1]);
      }
    } else {
      if (rank > 0) {
        //recv left
        MPI_Irecv(image, ny, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &requests[1]);
        //send left
        MPI_Isend(image + ny, ny, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &requests[0]);

        //recv right
        MPI_Irecv(image + ((width + 1) * ny), ny, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &requests[3]);
        //send right
        MPI_Isend(image + (width * ny), ny, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &requests[2]);
      } else {
        //recv right
        MPI_Irecv(image + ((width + 1) * ny), ny, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &requests[3]);
        //send right
        MPI_Isend(image + (width * ny), ny, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &requests[2]);
      }
    }

    // do the operations that don't rely on the completion of the nonblocking send-recvs first
    // to allow some time for the data transfer to complete
    for (int i = 2; i < width; ++i) {
      for (int j = 0; j < ny; ++j) {
        tmp_image[j+i*ny] = image[j+i*ny] * 0.6f;
        tmp_image[j + i * ny] += image[j + (i-1) * ny] * 0.1f;
        tmp_image[j + i * ny] += image[j + (i+1) * ny] * 0.1f;
        if (j > 0)          tmp_image[j+i*ny] += image[j-1+i*ny] * 0.1f;
        if (j < ny-1)       tmp_image[j+i*ny] += image[j+1+i*ny] * 0.1f;
      }
    }

    if (rank > 0) {
      MPI_Wait(&requests[1], MPI_STATUS_IGNORE);
      //MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
    }

    if (rank < size - 1) {
      MPI_Wait(&requests[3], MPI_STATUS_IGNORE);
      //MPI_Wait(&requests[2], MPI_STATUS_IGNORE);
    }

    int i = 1;
    for (int j = 0; j < ny; ++j) {
      tmp_image[j+i*ny] = image[j+i*ny] * 0.6f;
      if (i > 0)          tmp_image[j + i * ny] += image[j + (i-1) * ny] * 0.1f;
      tmp_image[j+i*ny] += image[j  +(i+1)*ny] * 0.1f;
      if (j > 0)          tmp_image[j+i*ny] += image[j-1+i*ny] * 0.1f;
      if (j < ny-1)       tmp_image[j+i*ny] += image[j+1+i*ny] * 0.1f;
    }

    i = width;
    for (int j = 0; j < ny; ++j) {
      tmp_image[j+i*ny] = image[j+i*ny] * 0.6f;
      tmp_image[j + i * ny] += image[j + (i-1) * ny] * 0.1f;
      if (i < nx-1)       tmp_image[j+i*ny] += image[j  +(i+1)*ny] * 0.1f;
      if (j > 0)          tmp_image[j+i*ny] += image[j-1+i*ny] * 0.1f;
      if (j < ny-1)       tmp_image[j+i*ny] += image[j+1+i*ny] * 0.1f;
    }

}

void stencilSyncCols(const int nx, const int ny, float * restrict image, float * restrict tmp_image, int width) {

  if (rank > 0 && rank < size - 1) {
    // left
    MPI_Sendrecv(image + ny, ny, MPI_FLOAT, rank - 1, 0, image, ny, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //right
    MPI_Sendrecv(image + (width * ny), ny, MPI_FLOAT, rank + 1, 0, image + ((width + 1) * ny), ny, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else if (rank > 0) {
    // left
    MPI_Sendrecv(image + ny, ny, MPI_FLOAT, rank - 1, 0, image, ny, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    //right
    MPI_Sendrecv(image + (width * ny), ny, MPI_FLOAT, rank + 1, 0, image + ((width + 1) * ny), ny, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  for (int i = 1; i <= width; ++i) {
    for (int j = 0; j < ny; ++j) {
      tmp_image[j+i*ny] = image[j+i*ny] * 0.6f;
      if (i > 0)          tmp_image[j+i*ny] += image[j  +(i-1)*ny] * 0.1f;
      if (i < nx-1)       tmp_image[j+i*ny] += image[j  +(i+1)*ny] * 0.1f;
      if (j > 0)          tmp_image[j+i*ny] += image[j-1+i*ny] * 0.1f;
      if (j < ny-1)       tmp_image[j+i*ny] += image[j+1+i*ny] * 0.1f;
    }
  }
}

// Create the input image
void init_image(const int nx, const int ny, float *  image, float *  tmp_image) {
  // Zero everything
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      image[j+i*ny] = 0.0;
      tmp_image[j+i*ny] = 0.0;
    }
  }

  // Checkerboard
  for (int j = 0; j < 8; ++j) {
    for (int i = 0; i < 8; ++i) {
      for (int jj = j*ny/8; jj < (j+1)*ny/8; ++jj) {
        for (int ii = i*nx/8; ii < (i+1)*nx/8; ++ii) {
          if ((i+j)%2)
          image[jj+ii*ny] = 100.0;
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char * file_name, const int nx, const int ny, float *image) {

  // Open output file
  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  float maximum = 0.0;
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (image[j+i*ny] > maximum)
        maximum = image[j+i*ny];
    }
  }

  // Output image, converting to numbers 0-255
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      fputc((char)(255.0*image[j+i*ny]/maximum), fp);
    }
  }

  // Close the file
  fclose(fp);

}

// Get the current time in seconds since the Epoch
double  wtime(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}
