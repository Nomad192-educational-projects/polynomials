#ifndef POLY_H
#define POLY_H

#include <stdio.h>

#define POLY_TYPE double

typedef struct poly_t
{
    size_t n;
    POLY_TYPE *data;
} Poly;

void print_poly(Poly the_poly);

Poly add_poly(Poly first, Poly second);
Poly sub_poly(Poly first, Poly second);
Poly mul_poly (Poly first, Poly second);

typedef struct res_dvn_t
{
    Poly remains;
    Poly result;
} Res_Dvn;

Res_Dvn dvn_poly (Poly divisible, Poly divider);
#endif