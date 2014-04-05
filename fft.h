#ifndef FFT_H
#define FFT_H
#include<math.h>
const float pi = 3.141592653589793238462643;
class vec2{
public:
    float x, y;
    vec2(float _x, float _y){
        x = _x, y = _y;
    }
    vec2(float _x){
        x = _x, y = 0;
    }
    vec2(){
        x = 0, y = 0;
    }
    vec2 operator *(vec2 r){
        return vec2(x * r.x - y * r.y, x * r.y + y * r.x);
    }
    vec2 operator *(float r){
        return vec2(x * r, y * r);
    }
    vec2 operator +(vec2 r){
        return vec2(x + r.x, y + r.y);
    }
    vec2 operator -(vec2 r){
        return vec2(x - r.x, y - r.y);
    }
    float r2(){
        return x * x + y * y;
    }
};
float abs(vec2 a){
    return sqrt(a.x*a.x+a.y*a.y);
}
void fft(float *a, vec2 *y, int n, int step = 1){
    //static vec2 yb[1024*1024];
    if(n == 1){
        *y = *a;
        return;
    }
    int k = n >> 1;
    fft(a, y, k, step << 1);
    fft(a + step, y + k, k, step << 1);
    vec2 w(1),  wn(cosf(2*pi/n), sinf(2*pi/n));
    for (int i = 0; i < k; ++i) {
        vec2 yi = y[i] + w * y[i+k], yik = y[i] - w * y[i+k];
        y[i] = yi;
        y[i+k] = yik;
        w = w * wn;
    }
}
void fft(vec2 *a, vec2 *y, int n, int step = 1){
    //static vec2 yb[1024*1024];
    if(n == 1){
        *y = *a;
        return;
    }
    int k = n >> 1;
    fft(a, y, k, step << 1);
    fft(a + step, y + k, k, step << 1);
    vec2 w(1),  wn(cosf(2*pi/n), sinf(2*pi/n));
    for (int i = 0; i < k; ++i) {
        vec2 yi = y[i] + w * y[i+k], yik = y[i] - w * y[i+k];
        y[i] = yi;
        y[i+k] = yik;
        w = w * wn;
    }
}
void ifft(vec2 *a, vec2 *y, int n, int step = 1){
    //static vec2 yb[1024*1024];
    if(n == 1){
        *y = *a;
        return;
    }
    int k = n >> 1;
    ifft(a, y, k, step << 1);
    ifft(a + step, y + k, k, step << 1);
    vec2 w(1),  wn(cosf(2*pi/n), -sinf(2*pi/n));
    for (int i = 0; i < k; ++i) {
        vec2 yi = y[i] + w * y[i+k], yik = y[i] - w * y[i+k];
        y[i] = yi*0.5f;
        y[i+k] = yik*0.5f;
        w = w * wn;
    }
}
#endif // FFT_H
