#ifndef __sysutil__
#define __sysutil__

int add(float* src, float* dest)
{
    // puts result in dest
    for(int i=0; i< ndim; i++)
        dest[i] = dest[i] + src[i];
    return 0;
}
int copy(float* src, float* dest)
{
    for(int i=0; i< ndim; i++)
        dest[i] = src[i];
    return 0;
}
float max_norm(float* src, int dim=ndim)
{
    float m = -1;
    for(int i=0; i<dim; i++)
    {
        if(m < src[i])
            m = src[i];
    }
    return m;
}
int print_state(float* src, int dim=ndim)
{
    for(int i=0; i<dim; i++)
        cout<<src[i]<<" ";
    cout<<endl;
    return 0;
}
float norm(float* s)
{
    float sum = 0;
    for(int i=0; i<ndim;i++)
        sum = sum + sq(s[i]);
    return sqrt(sum);
}
#endif
