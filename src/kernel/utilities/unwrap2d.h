#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// TODO: remove global variables
// TODO: make thresholds independent

static float PI = 3.141592654f;
static float TWOPI = 6.283185307f;

#define NOMASK 0
#define MASK 1

typedef struct {
  float mod;
  int x_connectivity;
  int y_connectivity;
  int no_of_edges;
} params_t;

// PIXELM information
struct PIXELM {
  int increment;  // No. of 2*pi to add to the pixel to unwrap it
  int number_of_pixels_in_group;  // No. of pixel in the pixel group
  float value;                    // value of the pixel
  float reliability;
  unsigned char input_mask;     // 0 pixel is masked. NOMASK pixel is not masked
  unsigned char extended_mask;  // 0 pixel is masked. NOMASK pixel is not masked
  int group;                    // group No.
  int new_group;
  struct PIXELM
      *head;  // pointer to the first pixel in the group in the linked list
  struct PIXELM *last;  // pointer to the last pixel in the group
  struct PIXELM *next;  // pointer to the next pixel in the group
};

typedef struct PIXELM PIXELM;

// the EDGE is the line that connects two pixels.
// if we have S pixels, then we have S horizontal edges and S vertical edges
struct EDGE {
  float reliab;       // reliabilty of the edge and it depends on the two pixels
  PIXELM *pointer_1;  // pointer to the first pixel
  PIXELM *pointer_2;  // pointer to the second pixel
  int increment;      // No. of 2*pi to add to one of the pixels to
                  // unwrap it with respect to the second
};

typedef struct EDGE EDGE;

//---------------start quicker_sort algorithm --------------------------------
#define swap(x, y) \
  {                \
    EDGE t;        \
    t = x;         \
    x = y;         \
    y = t;         \
  }
#define order(x, y) \
  if (x.reliab > y.reliab) swap(x, y)
#define o2(x, y) order(x, y)
#define o3(x, y, z) \
  o2(x, y);         \
  o2(x, z);         \
  o2(y, z)

typedef enum { yes, no } yes_no;




yes_no find_pivot(EDGE *left, EDGE *right, double *pivot_ptr);

EDGE *partition(EDGE *left, EDGE *right, double pivot);

void quicker_sort(EDGE *left, EDGE *right);
//--------------end quicker_sort algorithm -----------------------------------

//--------------------start initialize pixels ----------------------------------
// initialize pixels. See the explination of the pixel class above.
// initially every pixel is assumed to belong to a group consisting of only
// itself
void initialisePIXELs(double *wrapped_image, unsigned char *input_mask,
                      unsigned char *extended_mask, PIXELM *pixel,
                      int image_width, int image_height);
//-------------------end initialize pixels -----------

// gamma function in the paper
double wrap(double pixel_value);

// pixelL_value is the left pixel,	pixelR_value is the right pixel
int find_wrap(double pixelL_value, double pixelR_value);

void extend_mask(unsigned char *input_mask, unsigned char *extended_mask,
                 int image_width, int image_height, params_t *params);

void calculate_reliability(double *wrappedImage, PIXELM *pixel, int image_width,
                           int image_height, params_t *params);

// calculate the reliability of the horizontal edges of the image
// it is calculated by adding the reliability of pixel and the relibility of
// its right-hand neighbour
// edge is calculated between a pixel and its next neighbour
void horizontalEDGEs(PIXELM *pixel, EDGE *edge, int image_width,
                     int image_height, params_t *params);
// calculate the reliability of the vertical edges of the image
// it is calculated by adding the reliability of pixel and the relibility of
// its lower neighbour in the image.
void verticalEDGEs(PIXELM *pixel, EDGE *edge, int image_width, int image_height,
                   params_t *params);

// gather the pixels of the image into groups
void gatherPIXELs(EDGE *edge, params_t *params);

// unwrap the image
void unwrapImage(PIXELM *pixel, int image_width, int image_height);
// set the masked pixels (mask = 0) to the minimum of the unwrapper phase
void maskImage(PIXELM *pixel, unsigned char *input_mask, int image_width,
               int image_height);

// the input to this unwrapper is an array that contains the wrapped
// phase map.  copy the image on the buffer passed to this unwrapper to
// over-write the unwrapped phase map on the buffer of the wrapped
// phase map.
void returnImage(PIXELM *pixel, double *unwrapped_image, int image_width,
                 int image_height);

// the main function of the unwrapper
void unwrap2D(float *wrapped_image, float *UnwrappedImage,
              unsigned char *input_mask, int image_width, int image_height,
              int wrap_around_x, int wrap_around_y);
