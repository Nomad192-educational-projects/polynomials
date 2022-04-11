#include "poly.h"

static POLY_TYPE my_abs(POLY_TYPE a)
{
    if (a < 0)
        a *= -1.0;
    return a;
}

static void print_monomial(const char *sign, POLY_TYPE val, size_t i)
{
    if(i > 1)
        printf("%s%"SP"x^%zu",  sign, my_abs(val), i);
    else if(i == 1)
        printf("%s%"SP"x",      sign, my_abs(val));
    else
        printf("%s%"SP"",       sign, my_abs(val));
}

void print_poly(Poly the_poly)
{
    size_t i = the_poly.n;
    while (i--)
        if(my_abs(the_poly.data[i]) > 6E-7 || the_poly.n == 1 || i == 0)
        {
            print_monomial((the_poly.data[i] >= 0.0 ? "" : "- "), the_poly.data[i],  i);
            break;
        }
    while (i--)
        if(my_abs(the_poly.data[i]) > 6E-7)
            print_monomial((the_poly.data[i] >= 0.0 ? " + " : " - "), the_poly.data[i], i);
    printf("\n");
}

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))

static POLY_TYPE add(POLY_TYPE first, POLY_TYPE second)
{
    return first + second;
}

static POLY_TYPE sub(POLY_TYPE first, POLY_TYPE second)
{
    return first - second;
}

static Poly op_poly(POLY_TYPE (*f)(POLY_TYPE, POLY_TYPE), Poly first, Poly second)
{
    const size_t poly_res_len = MAX(first.n, second.n);
    Poly result = {poly_res_len, malloc(poly_res_len * sizeof(POLY_TYPE))};
    if(result.data == NULL)
        return result;

    size_t i = 0;
    for(; i < MIN(first.n, second.n); i++)
    {
        result.data[i] = f(first.data[i], second.data[i]);
    }
    if (first.n > second.n)
        for(; i < poly_res_len; i++)
            result.data[i] = first.data[i];
    else
        for(; i < poly_res_len; i++)
            result.data[i] =  second.data[i];

    return result;
}

Poly add_poly(Poly first, Poly second)
{
    return op_poly(add, first, second);
}

Poly sub_poly(Poly first, Poly second)
{
    return op_poly(sub, first, second);
}

Poly mul_poly (Poly first, Poly second)
{
    const size_t poly_res_len = (first.n + second.n) - 1;
    Poly result = {poly_res_len, malloc(poly_res_len * sizeof(POLY_TYPE))};
    if(result.data == NULL)
        return result;

    for(size_t i = 0; i < poly_res_len; i++)
        result.data[i] = 0;

    for(size_t i = 0; i < first.n; i++)
        for(size_t j = 0; j < second.n; j++)
            result.data[i + j] += first.data[i] * second.data[j];

    return result;
}

static _Bool Eq(double a, double b)
{
    return (my_abs(a - b) <= MAX(my_abs(a), my_abs(b)) * 6E-7);
}

Res_Dvn dvn_poly (Poly divisible, Poly divider)
{
    Res_Dvn res = {{0, NULL}, {0, NULL}};
    Poly result;
    Poly remains;
    {
        size_t len_res = (divisible.n - divider.n + 1);
        size_t len_rem = divisible.n;
        Poly buf_res = {len_res, (POLY_TYPE *)malloc(len_res * sizeof(POLY_TYPE))};
        if(buf_res.data == NULL)
            return res;
        Poly buf_rem = {len_rem, (POLY_TYPE *)malloc(len_rem * sizeof(POLY_TYPE))};
        if(buf_rem.data == NULL)
        {
            free(buf_res.data);
            return res;
        }
        result = buf_res;
        remains = buf_rem;

        for(size_t i = 0; i < remains.n; i++)
            remains.data[i] = divisible.data[i];

        for(size_t i = 0; i < result.n; i++)
            result.data[i] = 0;
    }



    size_t n = 0;

    while(1)
    {
        POLY_TYPE ratio = remains.data[remains.n - n - 1] / divider.data[divider.n - 1];

        remains.data[remains.n - n - 1] = 0.0;

        for(size_t i = remains.n - n - 1, j = divider.n - 1; i-- && j--;)
        {
            if(Eq(remains.data[i], divider.data[j] * ratio))
                n++;
            remains.data[i] -= divider.data[j] * ratio;
        }
        result.data[remains.n - n - divider.n] = ratio;
        n++;
        //print_poly(res_dvn.remains);
        if(divisible.n - n < divider.n)
        {
            POLY_TYPE *new_data = realloc(remains.data, (remains.n - n) * sizeof(POLY_TYPE));
            if (new_data == NULL)
            {
                free(result.data);
                free(remains.data);
                return res;
            }
            remains.data = new_data;
            remains.n -= n;
            //print_poly(res_dvn.remains);

            res.result = result;
            res.remains = remains;
            return res;
        }
    }
}