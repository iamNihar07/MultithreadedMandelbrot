
#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <time.h>

struct parameter //structure to store essential data per each thread
{
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	//double scale;
	//int num;
	int max;
	int thread_id;
};

int iteration_to_color(int i, int max);
int iterations_at_point(double x, double y, int max);
void *compute_image(void *ptr); //thread function which will be used
void compute_image1(struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max); //function used when n is not input or n=1. also serves as a good way to debug
struct bitmap *bm; //defining global so that i dont have to pass it everytime to functions
int num; //stores the number of threads


void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-n <thread> Specify the number of threads you want to be running\n"); //added this printf statement
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000 -n 50\n\n");
}

int main(int argc, char *argv[])
{	

	struct timespec start, finish;
	double elapsed;

	clock_gettime(CLOCK_MONOTONIC, &start); //getting the current time value for beginning 


	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int image_width = 500;
	int image_height = 500;
	int max = 1000;
	int n = 1; //Stores the number of threads and is 1 by default.

	// For each command line argument given,
	// override the appropriate configuration value.

	while ((c = getopt(argc, argv, "x:y:s:W:H:m:o:n:h")) != -1) //added command line n here
	{
		switch (c)
		{
		case 'x':
			xcenter = atof(optarg);
			break;
		case 'y':
			ycenter = atof(optarg);
			break;
		case 's':
			scale = atof(optarg);
			break;
		case 'W':
			image_width = atoi(optarg);
			break;
		case 'H':
			image_height = atoi(optarg);
			break;
		case 'm':
			max = atoi(optarg);
			break;
		case 'o':
			outfile = optarg;
			break;
		case 'n': // in case arg is -n
			n = atoi(optarg); //convert the string argument to number and store it in n(number of threads)
			//printf("%d\n", n);
			break;
		case 'h':
			show_help();
			exit(1);
			break;
		}
	}

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s threads=%d\n", xcenter, ycenter, scale, max, outfile, n);

	struct parameter instance[n]; //create n stances of the above defined structure for n threads

	// Create a bitmap of the appropriate size.
	bm = bitmap_create(image_width, image_height);

	num = n; //setting global var = local var so that i dont have to pass it to all functions

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm, MAKE_RGBA(0, 0, 255, 0));

	// Compute the Mandelbrot image
	if (n == 1) //serves as a debugging process and works for one process or thread =1 
	{
		compute_image1(bm, xcenter - scale, xcenter + scale, ycenter - scale, ycenter + scale, max);
	}
	else // if thread is more than one (multithreaded)
	{
	

		int i = 0; //counter variable for number of threads
		pthread_t thread[n]; //define n threads
		for (i = 0; i < n; i++)
		{
			//initializing each thread's attributes using the above defined struct
			instance[i].xmin = xcenter - scale;
			instance[i].xmax = xcenter + scale;
			instance[i].ymin = ycenter - scale;
			instance[i].ymax = ycenter + scale;
			instance[i].max = max;
			instance[i].thread_id = i; //storing the thread id (unique)
			if (pthread_create(&thread[i], NULL, compute_image, (void *)&instance[i])) //create n threads in for loop
			{
				perror("Exit creating threads");
				exit(EXIT_FAILURE);
			}			
		}

		for (i = 0; i < n; i++)
		{
			if (pthread_join(thread[i], NULL)) //synchronization of the created n threads
			{
				perror("Exit joining threads");
				exit(EXIT_FAILURE);
			}
		}
		
	}

	// Save the image in the stated file.
	if (!bitmap_save(bm, outfile))
	{
		fprintf(stderr, "mandel: couldn't write to %s: %s\n", outfile, strerror(errno));
		return 1;
	}												

	clock_gettime(CLOCK_MONOTONIC, &finish); //getting the time after threads have drawn the bitmap
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0; //getting the time elapsed from nano secs to secs.

	printf("Execution Time %0.02f sec.\n", elapsed); //printing elapsed time to check efficiency of multithreads.

	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void *compute_image(void *temp) //multithreading function
{

	//re obtain the values from the formal parameter passed on to this function

	struct parameter *ptr = (struct parameter *)temp;
	//printf("hi %d\n", ptr->thread_id);

	double xmin = ptr->xmin;
	double xmax = ptr->xmax;
	double ymin = ptr->ymin;
	double ymax = ptr->ymax;
	int max = ptr->max;
	int t_id = ptr->thread_id; //thread id 

	int i, j;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

	int divisions = height / num; //find out how many divisions of the pixels height we need in order to assign each thread the same (or almost) same amount of work
	int height1 = t_id * divisions; //serves as a lower bound for the y-limit when drawing the image 
	int height2 = t_id * divisions + divisions; //serves as a upper bound for the y-limit when drawing the image
	if(height2 >= height) //if the last allocation of threads exceeds the number of pixels, set it back to the max number of pixels.
	{					// this can happen when number of pixels are not divisibly by number of threads
		height2 = height-1;
	}
	//printf("%d %d %d %d \n",height, height1,height2, t_id);
	

	// For every pixel in the image...
	for (j = height1; j < height2; j++) //setting our new lower and upper bounds for each thread according to their id
	{

		//printf("%d \n", j);
		for (i = 0; i < width; i++)
		{

			// Determine the point in x,y space for that pixel.
			double x = xmin + i * (xmax - xmin) / width;
			double y = ymin + j * (ymax - ymin) / height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x, y, max);

			// Set the pixel in the bitmap.
			bitmap_set(bm, i, j, iters);
		}
	}
	
	//printf("bye %d\n", ptr->thread_id);
	
	pthread_exit(NULL); //exiting each thread when its job is done
}

void compute_image1(struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max) //as is, copied from the Github CSE 3320 repo, works for n=1 or unspecified n
{
	int i, j;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

	// For every pixel in the image...

	for (j = 0; j < height; j++)
	{

		for (i = 0; i < width; i++)
		{

			// Determine the point in x,y space for that pixel.
			double x = xmin + i * (xmax - xmin) / width;
			double y = ymin + j * (ymax - ymin) / height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x, y, max);

			// Set the pixel in the bitmap.
			bitmap_set(bm, i, j, iters);
		}
	}
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point(double x, double y, int max) //unchanged from cse 3320 github repo
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while ((x * x + y * y <= 4) && iter < max)
	{

		double xt = x * x - y * y + x0;
		double yt = 2 * x * y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter, max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color(int i, int max) //unchanged
{
	int gray = 255 * i / max;
	return MAKE_RGBA(gray, gray, gray, 0);
}
