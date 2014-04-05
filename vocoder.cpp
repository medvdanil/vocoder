#include<stdio.h>
#include<stdlib.h>
#include<memory.h>
#include"fft.h"
int stftsz = 512;
int smooth = 1;
const int wstep = stftsz/4;
struct IDRiff{
    char id[4];
    long len;
};


struct IDChuckWave{
    char id[4];
    char fmt[4];
    long len;
};

struct IDWave{
        short type;
        short channels;
        long SamplesPerSec;
        long AvgBytesPerSec;
        short align;
        short bits;
        short cbsize;
        char ext[46];
};
struct IDSampleWave{
    char id[4];
    long sz4;
    long len;
};

struct IDData{
    char id[4];
    long len;
};
struct smpdata{
    long n;
    float *cl, *cr;
    smpdata(int N = 0){
        n = N;
        cl = cr = 0;
    }
    ~smpdata(){
        if(cl)
            free(cl);
        if(cr && cr != cl)
            free(cr);
    }

};
int *data[2];
char buff[256];
bool wavopen(const char* fname, smpdata &d){
    FILE *f = fopen(fname, "rb+");
    IDRiff riff;
    IDChuckWave chunk;
    IDWave w;
    if(fread(&riff, sizeof(IDRiff), 1, f) <= 0 ||
            !(riff.id[0] == 'R' && riff.id[1] == 'I' && riff.id[2] == 'F' && riff.id[3] == 'F')){
        printf("It is not a WAV file\n");
        return 0;
    }
    fread(&chunk, sizeof(IDChuckWave), 1, f);
    fread(&w, chunk.len, 1, f);

    IDSampleWave smp;
    IDData dt;
    fread(buff, 4, 1, f);
    if(buff[0] == 'f'){
        fread(smp.id+4, sizeof(IDSampleWave)-4, 1, f);
        fread(&dt, sizeof(IDData), 1, f);
    }
    else
        fread(dt.id+4, sizeof(IDData)-4, 1, f);
    d.n = dt.len / w.align;
    if(d.n <= 0)
        return 0;
    long smpsz = w.align/w.channels;
    d.cl = new float[d.n];
    if(w.channels == 2)
        d.cr = new float[d.n];
    else d.cr = 0;
    unsigned char *bf = new unsigned char[dt.len];
    fread(bf, dt.len, 1, f);
    if(w.channels == 2){
        for(int i = 0; i < d.n; i++){
            d.cl[i] = 0,   d.cr[i] = 0;
            if(smpsz == 1)
                d.cl[i] = ((float)bf[i*w.channels]-128)/128.0f,
                    d.cr[i] = ((float)bf[i*w.channels+1]-128)/128.0f;
            if(smpsz == 2)
                d.cl[i] = ((short*)bf)[i*w.channels]/32768.0f, d.cr[i] = ((short*)bf)[i*w.channels+1]/32768.0f;
        }
    }
    else{
        for(int i = 0; i < d.n; i++){
            d.cl[i] = 0;
            if(smpsz == 1)
                d.cl[i] = ((float)bf[i*w.channels]-128)/128.0f;
            if(smpsz == 2)
                d.cl[i] = ((short*)bf)[i*w.channels]/32768.0f;
        }
    }
    delete[] bf;
    fclose(f);
    return 1;
}
bool wavout(const char* fname, smpdata &d, int freq){
    FILE *f = fopen(fname, "wb+");
    short nc = 2 - (d.cr == 0);
    short typesz = sizeof(short);
    IDRiff riff = {{'R', 'I', 'F', 'F'}, (long)sizeof(IDChuckWave)+16+(long)sizeof(IDData) + d.n*nc*typesz};
    IDChuckWave chunk = {{'W', 'A', 'V', 'E'}, {'f', 'm', 't', ' '}, 16};
    IDWave w = {0x1, nc, freq, freq*nc*typesz, short(nc*typesz), short(8*typesz), 0, ""};
    IDData dt = { {'d', 'a', 't', 'a'}, nc*d.n*typesz};
    fwrite(&riff, sizeof(IDRiff), 1, f);
    fwrite(&chunk, sizeof(IDChuckWave), 1, f);
    fwrite(&w, chunk.len, 1, f);
    fwrite(&dt, sizeof(IDData), 1, f);
    short *bf = new short[d.n*nc];
    if(nc == 1)
        for(int i = 0; i < d.n; i++)
            bf[i] = d.cl[i] < -1 ? -32768 : (d.cl[i] < 1 ? d.cl[i]*32768.0f : 32767);
    else
        for(int i = 0; i < d.n; i++){
            bf[i*2] = d.cl[i] < -1 ? -32768 : (d.cl[i] < 1 ? d.cl[i]*32768.0f : 32767);
            bf[i*2+1] = d.cr[i] < -1 ? -32768 : (d.cr[i] < 1 ? d.cr[i]*32768.0f : 32767);
        }
    fwrite(bf, typesz*nc, d.n, f);
    fclose(f);
    return 1;
}
void swap(int &a, int &b){
    int c = a;
    a = b, b = c;
}
void swap(float &a, float &b){
    float c = a;
    a = b, b = c;
}
void gen(smpdata &d, float *a, int n){
    d.cl = new float[d.n];
    d.cr = 0;
    srand(d.cl-d.cr);
    for(int i = 0; i < d.n; i++){
        d.cl[i] = 0;
        for(int j = 0; j < n; j++)
            d.cl[i] += sin(i*M_PI*a[j]/44100.0f);
        d.cl[i]/=n;
    }
}
int I = 0;
float* merge(float* a, float* b, int n){
    //    static float wa[stftsz*2];
     //   float wb[stftsz*2];
      //  static vec2 ya[stftsz*2];
       // static vec2 yb[stftsz*2];
    float wa[stftsz*2];
    float wb[stftsz*2];
    vec2 ya[stftsz*2];
    vec2 yb[stftsz*2];
    vec2 yh[stftsz*2];
    if(!a || !b) return 0;
    float *r = (float*)malloc(n*sizeof(float));
    memset(wa, 0 , sizeof(float)*stftsz);
    memset(wb, 0 , sizeof(float)*stftsz);
    memset(r, 0 , sizeof(float)*n);
    for(int i = -stftsz+wstep; i < n; i += wstep){
        if(i*80/n > (i-wstep)*80/n)
            putchar('=');
        for(int j = 0; j < stftsz; j++)
            if(i+j < 0 || i+j >= n)
                wa[j] = 0, wb[j] = 0;
            else{
                float h  = (1 - cos(2 * pi * (j + 0.5f) / stftsz))/(stftsz/wstep);
                wa[j] = h * a[i+j],
                    wb[j] = h * b[i+j];
            }

        fft(wa, ya, stftsz);
        fft(wb, yb, stftsz);

        //for(int j = 0; j < stftsz; j++)
           //ya[j] = ((yb[j].r2() < eps) ? vec2(0, 0) : yb[j] * sqrt(ya[j].r2() / yb[j].r2()));

        for(int j = 0; j < stftsz; j++)
            wb[j] = ((yb[j].r2()  == 0) ? 1 : sqrt(ya[j].r2() / yb[j].r2()));

       for(int j = 0; j < stftsz; j++){
            wa[j] = wb[j];
            int r = smooth;
            double s = 1;
            for(int k = 1; k < r; k++){
                double h = (1 - sin(0.5*pi * k / r));
                s+=h*2;
                wa[j] += h*( (j-k >= 0 ? wb[j-k] : 0)
                        + (j+k < stftsz ? wb[j+k] : 0));
            }
            wa[j]/=s;
        }
        for(int j = 0; j < stftsz; j++)
            ya[j] = yb[j] * wa[j];

        ifft(ya, yb, stftsz);
        for(int j = 0; j < stftsz; j++)
            if(i+j >= 0 && i+j < n)
                r[i+j] += yb[j].x*(1 - cos(2 * pi * (j + 0.5f) / stftsz));
    }
    return r;
}
float* extend(float* a, int n){
    float wa[stftsz];
    vec2 ya[stftsz];
    vec2 yr[stftsz];
    if(!a) return 0;
    float *r = (float*)malloc(n*sizeof(float)*2);
    memset(wa, 0 , sizeof(float)*stftsz);
    memset(r, 0, sizeof(float)*n*2);
    for(int i = -stftsz+wstep; i < n; i += wstep){
        if(i*80/n > (i-wstep)*80/n)
            putchar('=');
        for(int j = 0; j < stftsz; j++)
            if(i+j < 0 || i+j >= n)
                wa[j] = 0;
            else{
                float h  = (1 - cos(2 * pi * (j + 0.5f) / stftsz))/(stftsz/wstep);
                wa[j] = h * a[i+j];
            }
        for(int j = 0; j < stftsz/2; j++)
            swap(wa[j], wa[j+stftsz/2]);


        fft(wa, ya, stftsz);

        for(int j = 0; j < stftsz; j++){
                ya[j] = ((ya[j].r2() != 0) ? ya[j]*ya[j] * sqrt(1 / ya[j].r2()) : vec2(0, 0));
                      //  + ((ya[stftsz+j].r2() > eps) ? ya[stftsz+j]*ya[stftsz+j] * sqrt(1 / ya[stftsz+j].r2()) : vec2(0, 0));
        }
        ifft(ya, yr, stftsz);
        for(int j = 0; j < stftsz; j++)
            wa[j] = yr[j].x;//*(1 - cos(2 * pi * (j + 0.5f) / stftsz));
        for(int j = 0; j < stftsz/2; j++)
            swap(wa[j], wa[j+stftsz/2]);
        for(int j = 0; j < stftsz; j++)
            if(i*2+j >= 0 && i*2+j < n*2)
                r[i*2+j] += wa[j]*(1 - cos(2 * pi * (j + 0.5f) / stftsz));
    }
    return r;
}
float* freq(float* a, int n){
    float wa[stftsz];
    vec2 ya[stftsz];
    vec2 yr[stftsz];
    if(!a) return 0;
    float *r = (float*)malloc(n*sizeof(float));
    memset(wa, 0 , sizeof(float)*stftsz);
    memset(r, 0, sizeof(float)*n);
    for(int i = -stftsz+wstep; i < n; i += wstep){
        if(i*80/n > (i-wstep)*80/n)
            putchar('=');
        for(int j = 0; j < stftsz; j++)
            if(i+j < 0 || i+j >= n)
                wa[j] = 0;
            else{
                float h  = (1 - cos(2 * pi * (j + 0.5f) / stftsz))/(stftsz/wstep);
                wa[j] = h * a[i+j];
            }
        for(int j = 0; j < stftsz/2; j++)
            swap(wa[j], wa[j+stftsz/2]);


        fft(wa, ya, stftsz);

        for(int j = 0; j < stftsz; j++){
                ya[j] = ((ya[j].r2() != 0) ? ya[j]*ya[j] * sqrt(1 / ya[j].r2()) : vec2(0, 0));
                      //  + ((ya[stftsz+j].r2() > eps) ? ya[stftsz+j]*ya[stftsz+j] * sqrt(1 / ya[stftsz+j].r2()) : vec2(0, 0));
        }
        ifft(ya, yr, stftsz);
        for(int j = 0; j < stftsz; j++)
            wa[j] = yr[j].x;//*(1 - cos(2 * pi * (j + 0.5f) / stftsz));
        for(int j = 0; j < stftsz/2; j++)
            swap(wa[j], wa[j+stftsz/2]);
        for(int j = 0; j < stftsz; j++)
            if(i+j >= 0 && i+j < n)
                r[i+j] += wa[j]*(1 - cos(2 * pi * (j + 0.5f) / stftsz));
    }
    return r;
}
float* minus(float* a, float* b, int n){
    int i;
    if(!a || !b) return 0;
    float *r = (float*)malloc(n*sizeof(float));
    for(i = 0; i < n; i++)
        r[i] = a[i]-b[i];
    return r;
}

int main(int argc, char ** argv){
    //freopen("out.txt", "w", stdout);

    smpdata d1, d2;
    smpdata r;
    if(argc > 1 &&  argv[1][0] == '-' && argv[1][1] == '0'){
        r.n = 44100*15;
        float *a = new float[argc-3];
        for(int i = 3; i < argc; i++)
            a[i-3] = atof(argv[i]);
        gen(r, a, argc-3);
        wavout(argc > 2 ? argv[2] : "sampleout.wav", r, 44100);

        return 0;
    }
    if(argc < 2){
        if(!wavopen("Speech.wav", d1) || !wavopen("Music.wav", d2))        return 0;
    }
    else if(argv[1][0] == '-' && argv[1][1] == '1'){
        if(!wavopen(argv[2], d1) || !wavopen(argv[3], d2))        return 0;
        if(argc > 6 && argv[5][0] == '-' && argv[5][1] == 'w'){
            stftsz = atoi(argv[6]);
            if(stftsz <= 0)
               stftsz = 512;
        }
        if(argc > 8 && argv[7][0] == '-' && argv[7][1] == 'r'){
            smooth = atoi(argv[8]);
            if(smooth <= 0)
                smooth = 1;
        }
       // if(!wavopen("Speech.wav", d1) || !wavopen("Music.wav", d2))        return 0;
    }
    else{
        if(!wavopen(argv[2], d1)) return 0;
        if(argc > 5 && argv[4][0] == '-' && argv[4][1] == 'w'){
            stftsz = atoi(argv[5]);
            if(stftsz <= 0)
                stftsz = 512;
        }
        if(argc > 7 && argv[6][0] == '-' && argv[6][1] == 'r'){
            smooth = atoi(argv[7]);
            if(smooth <= 0)
                smooth = 1;
        }
    }
        //if(!wavopen("Music.wav", d1))        return 0;
    int n = (d1.n + stftsz-1)/ stftsz * stftsz;
    d1.cl = (float*)realloc(d1.cl, n*sizeof(float));
    memset(d1.cl+d1.n, 0, (n - d1.n)*sizeof(float));
    if(d1.cr){
        d1.cr = (float*)realloc(d1.cr, n*sizeof(float));
        memset(d1.cr+d1.n, 0, (n - d1.n)*sizeof(float));
    }
    d1.n = n;
    if(argc > 2 && (argv[1][0] == '-' && argv[1][1] == '1')){
        r.n = d1.n;
        if(d2.n < d1.n){
            n = d2.n * (d1.n/d2.n + (d1.n%d2.n != 0));
            d2.cl = (float*)realloc(d2.cl, n*sizeof(float));
            if(d2.cr) d2.cr = (float*)realloc(d2.cr, n*sizeof(float));
            for(int i = d2.n; i < n; i+=d2.n){
                memcpy(d2.cl+i, d2.cl, sizeof(float)*d2.n);
                if(d2.cr)
                    memcpy(d2.cr+i, d2.cr, sizeof(float)*d2.n);
            }
            d2.n = d1.n;
        }
        r.cl = merge(d1.cl, d2.cl, d1.n);
        r.cr = merge(d1.cr, d2.cr, d1.n);
        wavout(argc > 4 ? argv[4] : "sampleout.wav", r, 44100);
    }
    else if(argc > 2 && (argv[1][0] == '-' && argv[1][1] == '2')){
        r.n = d1.n*2;
        r.cl = extend(d1.cl, d1.n);
        r.cr = extend(d1.cr, d1.n);
        wavout(argc > 3 ? argv[3] : "sampleout.wav", r, 44100);
    }
    else if(argc > 2 && (argv[1][0] == '-' && argv[1][1] == '3')){
        r.n = d1.n;
        r.cl = freq(d1.cl, d1.n);
        r.cr = freq(d1.cr, d1.n);
        wavout(argc > 3 ? argv[3] : "sampleout.wav", r, 44100);
    }
    else if(argc > 2 && (argv[1][0] == '-' && argv[1][1] == '4')){
        r.n = d1.n*2;
        r.cl = (float*)malloc(r.n*sizeof(float));
        if(d1.cr)
            r.cr = (float*)malloc(r.n*sizeof(float));
        for(int i = 0;i < r.n;i++){
            r.cl[i] = (i&1) ? (d1.cl[i/2]+d1.cl[i/2+(i!=r.n-1)])/2 : d1.cl[i/2];
            if(r.cr)
                    r.cr[i] = (i&1) ? (d1.cr[i/2]+d1.cr[i/2+(i!=r.n-1)])/2 : d1.cr[i/2];
        }
        wavout(argc > 3 ? argv[3] : "sampleout.wav", r, 44100);
    }
    else{
        r.n = d1.n;
        if(d2.n < d1.n){
            n = d2.n * (d1.n/d2.n + (d1.n%d2.n != 0));
            d2.cl = (float*)realloc(d2.cl, n*sizeof(float));
            if(d2.cr) d2.cr = (float*)realloc(d2.cr, n*sizeof(float));
            for(int i = d2.n; i < n; i+=d2.n){
                memcpy(d2.cl+i, d2.cl, sizeof(float)*d2.n);
                if(d2.cr)
                    memcpy(d2.cr+i, d2.cr, sizeof(float)*d2.n);
            }
            d2.n = d1.n;
        }
        r.cl = minus(d1.cl, d2.cl, d1.n);
        r.cr = minus(d1.cr, d2.cr, d1.n);
        wavout(argc > 4 ? argv[4] : "sampleout.wav", r, 44100);
    }
    return 0;
}
