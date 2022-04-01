#ifndef WAV2MP3_H
#define WAV2MP3_H

#define HUFFBITS unsigned long int
#define HTN	34
#define MXOFF	250

struct huffcodetab {
    unsigned int xlen; 
    unsigned int ylen;
    unsigned int linbits;
    unsigned int linmax;
    HUFFBITS *table;
    unsigned char *hlen;
};
extern struct huffcodetab ht[HTN];

extern int sfBandIndex[23];
extern int scfsi_band_long[5];
extern double enwindow[512];
extern struct {
    unsigned region0_count;
    unsigned region1_count;
} subdv_table[23]; 

extern double cos_l[18][36];
extern double filter[32][64];
extern double ca[8];
extern double cs[8];
extern double win[36];

void subband_initialise(); 
void mdct_initialise();

#endif
