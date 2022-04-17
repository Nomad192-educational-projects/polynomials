#include <stdio.h>
#include <stdlib.h>
#include "poly.h"

Poly ini_poly(size_t n, const POLY_TYPE ini_list[n])
{
    Poly poly = {n, (POLY_TYPE*)malloc(n * sizeof(POLY_TYPE))};
    if(poly.data == NULL)
        return poly;

    for (size_t i = 0; i < poly.n; i++)
        poly.data[poly.n - 1 - i] = ini_list[i];
    return poly;
}

int main() {
    Poly first;
    Poly second;

    int main_result = 1;
    {
        POLY_TYPE ini_list_1[] = {2, 0, -5, 0, -8, 0};
        first = ini_poly(sizeof(ini_list_1)/sizeof(ini_list_1[0]), ini_list_1);
        if (first.data == NULL)
            goto exit;
        POLY_TYPE ini_list_2[] = {0, 0, 0, 1, 3};
        second = ini_poly(sizeof(ini_list_2)/sizeof(ini_list_2[0]), ini_list_2);
        if (second.data == NULL)
            goto exit;
    }


    printf("Poly 1: ");
    print_poly(first);
    printf("Poly 2: ");
    print_poly(second);

    {
        Poly op_result = add_poly(first, second);
        if(op_result.data == NULL)
            goto exit;
        printf("add 1+2: ");
        print_poly(op_result);
        free(op_result.data);
    }
    {
        Poly op_result = sub_poly(first, second);
        if(op_result.data == NULL)
            goto exit;
        printf("sub 1-2: ");
        print_poly(op_result);
        free(op_result.data);
    }
    {
        Poly op_result = mul_poly(first, second);
        if(op_result.data == NULL)
            goto exit;
        printf("mul 1*2: ");
        print_poly(op_result);
        free(op_result.data);
    }
    {
        Res_Dvn res_dvn = dvn_poly(first, second);
        if(res_dvn.result.data == NULL || res_dvn.remains.data == NULL)
            goto exit;
        printf("dvn 1/2:\nres = ");
        print_poly(res_dvn.result);
        printf("rem = ");
        print_poly(res_dvn.remains);
        free(res_dvn.result.data);
        free(res_dvn.remains.data);
    }

    main_result = 0;
    exit:
    free(first.data);
    free(second.data);
    return main_result;
}
