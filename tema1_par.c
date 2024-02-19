// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

typedef struct {
    int id;
    int nr_threads;
    pthread_barrier_t *barrier;
    ppm_image *image;
    ppm_image **contour_map;
    ppm_image *scaled_image;
    unsigned char **grid;
} my_arg;

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Allocates memory for the grid
unsigned char **aloc_sample_grid(ppm_image *image) {
    int p = image->x / STEP;
    int q = image->y / STEP;

    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    return grid;
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
void sample_grid(ppm_image *image, unsigned char **grid ,int step_x, int step_y, unsigned char sigma, int id, int nr_threads) {
    int p = image->x / step_x;
    int q = image->y / step_y;

    // segment the size of the image based on the number of threads
    int start = id * p / nr_threads;
    int end = fmin((id + 1) * p / nr_threads, p);

    // parallelize outer loop
    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = image->data[i * step_x * image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }
    grid[p][q] = 0;

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = image->data[i * step_x * image->y + image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }

    // segment the size of the image based on the number of threads
    int start1 = id * q / nr_threads;
    int end1 = fmin((id + 1) * q / nr_threads, q);

    for (int j = start1; j < end1; j++) {
        ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(ppm_image *image, unsigned char **grid, ppm_image **contour_map, int step_x, int step_y, int id, int nr_threads) {
    int p = image->x / step_x;
    int q = image->y / step_y;

    int start = id * p / nr_threads;
    int end = fmin((id + 1) * p / nr_threads, p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image, contour_map[k], i * step_x, j * step_y);
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

// Allocates memory for the rescaled image if needed
ppm_image *aloc_rescaled_image(ppm_image *image) {
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        return image;
    }

    ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    new_image->x = RESCALE_X;
    new_image->y = RESCALE_Y;

    new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    return new_image;
}

// Rescales the image using the bicubic interpolation algorithm
void rescale_image(ppm_image *image, ppm_image *new_image, int id, int nr_threads) {
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        new_image = image;
        return;
    }

    uint8_t sample[3];

    // segment the size of the image based on the number of threads
    int start = id * new_image->x / nr_threads;
    int end = fmin((id + 1) * new_image->x / nr_threads, new_image->x);

    // parallelize outer loop
    for (int i = start; i < end; i++) {
        for (int j = 0; j < new_image->y; j++) {
            float u = (float)i / (float)(new_image->x - 1);
            float v = (float)j / (float)(new_image->y - 1);
            sample_bicubic(image, u, v, sample);

            new_image->data[i * new_image->y + j].red = sample[0];
            new_image->data[i * new_image->y + j].green = sample[1];
            new_image->data[i * new_image->y + j].blue = sample[2];
        }
    }
}

void *thread_function(void *args) {
    my_arg *data = (my_arg *)args;
    int id = data->id;
    int nr_threads = data->nr_threads;
    pthread_barrier_t *barrier = data->barrier;
    ppm_image *image = data->image;
    ppm_image **contour_map = data->contour_map;

    int step_x = STEP;
    int step_y = STEP;

    // 1. Rescale the image
    rescale_image(image, data->scaled_image, id, nr_threads);

    // wait for all threads to finish rescaling
    pthread_barrier_wait(barrier);

    // save the finished rescaled image in a variable
    ppm_image *scaled_image = data->scaled_image;

    // 2. Sample the grid
    sample_grid(scaled_image, data->grid, step_x, step_y, SIGMA, id, nr_threads);

    // wait for all threads to finish sampling
    pthread_barrier_wait(barrier);

    // save the finished grid in a variable
    unsigned char **grid = data->grid;

    // 3. March the squares
    march(scaled_image, grid, contour_map, step_x, step_y, id, nr_threads);

    pthread_exit(NULL);
}


int main(int argc, char *argv[]) {
    // declare the variables used
    int nr_threads;
    void *status;
    pthread_t *threads;
    my_arg *args;
    pthread_barrier_t barrier;
    int r;

    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    ppm_image *image = read_ppm(argv[1]);

    // 0. Initialize contour map
    ppm_image **contour_map = init_contour_map();

    // allocate memory for the scaled image
    ppm_image *scaled_image = aloc_rescaled_image(image);

    // allocate memory for the grid
    unsigned char **grid = aloc_sample_grid(scaled_image);

    nr_threads = atoi(argv[3]);

    // allocate memory for threads and args
    threads = (pthread_t *) malloc(nr_threads * sizeof(pthread_t));
    args = (my_arg *) malloc(nr_threads * sizeof(my_arg));

    // initialize barrier
    r = pthread_barrier_init(&barrier, NULL, nr_threads);

    if (r) {
        printf("Could not initialize barrier\n");
        exit(-1);
    }

    // create threads
    for (int i = 0; i < nr_threads; i++) {
        args[i].id = i;
        args[i].nr_threads = nr_threads;
        // all the threads share the same barrier, image, contour_map, scaled_image and grid
        args[i].barrier = &barrier;
        args[i].image = image;
        args[i].contour_map = contour_map;
        args[i].scaled_image = scaled_image;
        args[i].grid = grid;

        r = pthread_create(&threads[i], NULL, thread_function, &args[i]);

        if (r) {
            printf("Could not create thread %d\n", i);
            exit(-1);
        }

    }

    // wait for threads to finish
    for (int i = 0; i < nr_threads; i++) {
        r = pthread_join(threads[i], &status);

        if (r) {
            printf("Could not join thread %d\n", i);
            exit(-1);
        }
    }

    write_ppm(scaled_image, argv[2]);

    // in case the initial image was rescaled, it needs to be freed
    // if they point to the same address, it means the image was not rescaled
    if (scaled_image != image) {
        free(image->data);
        free(image);
    }

    // free resources
    free_resources(scaled_image, contour_map, grid, STEP);

    // destroy barrier
    r = pthread_barrier_destroy(&barrier);
    if (r) {
        printf("Could not destroy barrier\n");
        exit(-1);
    }

    free(threads);
    free(args);

    return 0;
}
