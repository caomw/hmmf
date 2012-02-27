#ifndef __sysutil__
#define __sysutil__

int add(double* src, double* dest)
{
    // puts result in dest
    for(int i=0; i< ndim; i++)
        dest[i] = dest[i] + src[i];
    return 0;
}
int copy(double* src, double* dest)
{
    for(int i=0; i< ndim; i++)
        dest[i] = src[i];
    return 0;
}
double max_norm(double* src, int dim=ndim)
{
    double m = -1;
    for(int i=0; i<dim; i++)
    {
        if(m < src[i])
            m = src[i];
    }
    return m;
}
int print(double* src, int dim=ndim)
{
    for(int i=0; i<dim; i++)
        cout<<src[i]<<" ";
    cout<<endl;
    return 0;
}
double norm(double* s)
{
    double sum = 0;
    for(int i=0; i<ndim;i++)
        sum = sum + sq(s[i]);
    return sqrt(sum);
}
#endif
