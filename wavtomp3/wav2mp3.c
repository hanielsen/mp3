#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "wav2mp3_tbl.h"

typedef struct {
    unsigned part2_3_length;
    unsigned big_values;
    unsigned count1;
    unsigned global_gain;
    unsigned scalefac_compress;
    unsigned window_switching_flag;
    unsigned table_select[3];
    unsigned region0_count;
    unsigned region1_count;
    unsigned preflag;
    unsigned scalefac_scale;
    unsigned count1table_select;
    unsigned part2_length;
    unsigned sfb_lmax;
    unsigned address1;
    unsigned address2;
    unsigned address3;
    double quantizerStepSize;
    unsigned slen[4];
} gr_info;

struct side_info_t {
    unsigned main_data_begin;
    unsigned private_bits;
    unsigned scfsi[4];
    gr_info gr[2];
} side_info;

FILE *wav_file;
FILE *mp3_file;

double now_que[512];
double win_que[512]; 

int pos = 0;

short buffer[1152];
//struct side_info_t side_info;
double sb_sample[3][18][32];
double mdct_freq[2][576];
int enc[2][576];

void bit_pos(int n) 
{
    pos = n;
}        
    
int bit_set(int n, unsigned x, unsigned char bitbuf[]) 
{
    unsigned b, c;
    for (int i=0; i<n; i++) 
    {
        b = (x>>(n-1-i)) & 1;
        c = 1 << (7-pos%8);
        if (b) 
        {
            bitbuf[pos/8] |= c;
        }
        else
        {
            bitbuf[pos/8] &= (255-c);
        }
        pos += 1;
    } 
    return pos;
}
 
int abs_and_sign(int *x) 
{
    if (*x > 0) 
    {
        return 0;
    }
    *x *= -1;
    return 1;
}

int ix_max(int *ix, unsigned int begin, unsigned int end) 
{
    int x;
    int max = 0;
    for (int i=begin; i<end; i++) 
    { 
        x = abs(ix[i]);
        if (x > max) 
        {
            max = x;
        }
    }
    return max;
}

int HuffmanCode(int table_select, int x, int y, unsigned int *code, 
                unsigned int *ext, int *cbits, int *xbits ) 
{
    unsigned signx, signy, linbitsx, linbitsy, linbits, xlen, ylen, idx;
    struct huffcodetab *h;

    *cbits = 0;
    *xbits = 0;
    *code  = 0;
    *ext   = 0;
    
    if (table_select == 0) return 0;
    
    signx = abs_and_sign( &x );
    signy = abs_and_sign( &y );
    h = &(ht[table_select]);
    xlen = h->xlen;
    ylen = h->ylen;
    linbits = h->linbits;
    linbitsx = linbitsy = 0;

    if (table_select > 15) 
    {
        if ( x > 14 ) 
        {
            linbitsx = x - 15;
            x = 15;
        }
        if (y > 14) 
        {
            linbitsy = y - 15;
            y = 15;
        }

        idx = x*ylen + y;
        *code  = h->table[idx];
        *cbits = h->hlen[idx];
        if (x > 14) 
        {
            *ext |= linbitsx;
            *xbits += linbits;
        }
        if (x != 0) 
        {
            *ext <<= 1;
            *ext |= signx;
            *xbits += 1;
        }
        if (y > 14) 
        {
            *ext <<= linbits;
            *ext |= linbitsy;
            *xbits += linbits;
        }
        if (y != 0) 
        {
            *ext <<= 1;
            *ext |= signy;
            *xbits += 1;
        }
    }
    else 
    { 
        idx = x*ylen + y;
        *code = h->table[idx];
        *cbits += h->hlen[idx];
        if (x != 0) 
        {
            *code <<= 1;
            *code |= signx;
            *cbits += 1;
        }
        if (y != 0) 
        {
            *code <<= 1;
            *code |= signy;
            *cbits += 1;
        }
    }
    return *cbits + *xbits;
}

int huffman_coder_count1(unsigned char* pph, struct huffcodetab *h, int v, int w, int x, int y) 
{
    HUFFBITS huffbits;
    unsigned int signv, signw, signx, signy, p;
    int len;
    int totalBits = 0;
    
    signv = abs_and_sign(&v);
    signw = abs_and_sign(&w);
    signx = abs_and_sign(&x);
    signy = abs_and_sign(&y);
    
    p = v + (w << 1) + (x << 2) + (y << 3);
    huffbits = h->table[p];
    len = h->hlen[p];
    bit_set(len, huffbits, pph);
    totalBits += len;
    if (v) 
    {
        bit_set(1, signv, pph);
        totalBits += 1;
    }
    if (w) 
    {
        bit_set(1, signw, pph);
        totalBits += 1;
    }
    if (x) 
    {
        bit_set(1, signx, pph);
        totalBits += 1;
    }
    if (y) 
    {
        bit_set(1, signy, pph);
        totalBits += 1;
    }
    return totalBits;
}

void Huffmancodebits(unsigned char pph[397], int ix[576], gr_info *gi) 
{
    int region1Start;
    int region2Start;
    int bigvalues, count1End;
    int v, w, x, y, bits, cbits, xbits, stuffingBits;
    unsigned int code, ext;
    struct huffcodetab *h;
    int r0, r1, r2, rt, *pr;
    int bitsWritten = 0;
    r0 = r1 = r2 = 0;
    bigvalues = gi->big_values <<1;
    unsigned scalefac_index = gi->region0_count + 1;
    region1Start = sfBandIndex[scalefac_index];
    scalefac_index += gi->region1_count + 1;
    region2Start = sfBandIndex[scalefac_index];
    for (int i=0; i<bigvalues; i += 2 ) 
    {
        unsigned tableindex = 100;
        if (i < region1Start) 
        {
            tableindex = gi->table_select[0];
            pr = &r0;
        }
        else if (i < region2Start) 
        {
            tableindex = gi->table_select[1];
            pr = &r1;
        }
        else 
        {
            tableindex = gi->table_select[2];
            pr = &r2;
        }
        h = &ht[tableindex];
        x = ix[i];
        y = ix[i + 1];
        if (tableindex) 
        {
            bits = HuffmanCode(tableindex, x, y, &code, &ext, &cbits, &xbits);
            bit_set(cbits, code, pph);
            bit_set(xbits, ext, pph);
            rt = bits;
            bitsWritten += rt;
            *pr += rt;
        }
        else 
        {
            *pr = 0;
        }
    }
 
    h = &ht[gi->count1table_select + 32];
    count1End = bigvalues + (gi->count1 <<2);
    for (int i = bigvalues; i < count1End; i += 4) 
    {
        v = ix[i];
        w = ix[i+1];
        x = ix[i+2];
        y = ix[i+3];
        bitsWritten += huffman_coder_count1(pph, h, v, w, x, y);
    }
    bit_pos(gi->part2_3_length);
}

void huff()
{
    unsigned char main_data[397];
    gr_info* cod_info;
    
    for (int i=0; i<397; i++) 
    {
        main_data[i]=0;
    }
    bit_pos(0);
    for(int gr=0; gr<2; gr++) 
    {   
        cod_info = (gr_info*) &(side_info.gr[gr]);     
        Huffmancodebits(main_data, enc[gr], cod_info);
    }
    fwrite(main_data, 1, 397, mp3_file); 
}

int count_bit(int *ix, unsigned int start, unsigned int end, unsigned int table) 
{
    unsigned linbits, ylen;
    int sum;
    int x, y;
    struct huffcodetab *h;

    if (!table)
    {
        return 0;
    }
    h   = &(ht[table]);
    sum = 0;

    ylen    = h->ylen;
    linbits = h->linbits;

    if (table>15) 
    {  
        for (int i=start; i<end; i+=2) 
        {
            x = ix[i];
            y = ix[i+1];
            if (x>14) 
            {
                x = 15;
                sum += linbits;
            }
            if (y>14) 
            {
                y = 15;
                sum += linbits;
            }

            sum += h->hlen[x*ylen + y];

            if (x)
            {
                sum++;
            }
            if (y)
            {
                sum++;
            }
        }
    }
    else 
    { 
        for (int i=start; i<end; i+=2) 
        {
            x = ix[i];
            y = ix[i+1];

            sum  += h->hlen[(x*ylen)+y];

            if (x!=0) 
            {
                sum++;
            }
            if (y!=0)
            { 
                sum++;
            }
        }
    }

    return sum;
}

int bigv_bitcount(int *ix, gr_info *gi) 
{
    int bits = 0;
    unsigned int table;
        
    if (table = gi->table_select[0])  /* region0 */ 
    {
        bits += count_bit(ix, 0, gi->address1, table);
    }
    if (table = gi->table_select[1])  /* region1 */ 
    {
        bits += count_bit(ix, gi->address1, gi->address2, table);
    }
    if (table = gi->table_select[2])  /* region2 */ 
    {
        bits += count_bit(ix, gi->address2, gi->address3, table);
    }
    return bits;
}

int new_choose_table(int *ix, unsigned int begin, unsigned int end) 
{
    int i, max;
    int choice[2];
    int sum[2];

    max = ix_max(ix, begin, end);
    if (!max)
    {
        return 0;
    }
    choice[0] = 0;
    choice[1] = 0;

    if (max<15) 
    {
        for (int i=0; i<=13; i++) 
        {
            if (ht[i].xlen > max) 
            {
                choice[0] = i;
                break;
            }
        }
        sum[0] = count_bit(ix, begin, end, choice[0]);

        switch (choice[0]) 
        {
            case 2:
                sum[1] = count_bit(ix, begin, end, 3);
                if (sum[1] <= sum[0])
                {
                    choice[0] = 3;
                }
                break;

            case 5:
                sum[1] = count_bit(ix, begin, end, 6);
                if (sum[1] <= sum[0])
                {
                    choice[0] = 6;
                }
                break;

            case 7:
                sum[1] = count_bit(ix, begin, end, 8);
                if (sum[1] <= sum[0]) 
                {
                    choice[0] = 8;
                    sum[0] = sum[1];
                }
                sum[1] = count_bit(ix, begin, end, 9);
                if (sum[1] <= sum[0])
                {
                    choice[0] = 9;
                }
                break;

            case 10:
                sum[1] = count_bit(ix, begin, end, 11);
                if (sum[1] <= sum[0]) 
                {
                    choice[0] = 11;
                    sum[0] = sum[1];
                }
                sum[1] = count_bit(ix, begin, end, 12);
                if (sum[1] <= sum[0])
                {
                    choice[0] = 12;
                }
                break;

            case 13:
                sum[1] = count_bit(ix, begin, end, 15);
                if (sum[1] <= sum[0])
                {
                    choice[0] = 15;
                }
                break;
        }
    }
    else 
    {
        max -= 15;
    
        for (int i=15; i<24; i++) 
        {
            if (ht[i].linmax>=max) 
            {
                choice[0] = i;
                break;
            }
        }
        for (int i=24; i<32; i++) 
        {
            if (ht[i].linmax>=max) 
            {
                choice[1] = i;
                break;
            }
        }
    
        sum[0] = count_bit(ix, begin, end, choice[0]);
        sum[1] = count_bit(ix, begin, end, choice[1]);
        if (sum[1] < sum[0])
        {
            choice[0] = choice[1];
        }
    }
    return choice[0];
}

void bigv_tab_select(int *ix, gr_info *cod_info) 
{
    cod_info->table_select[0] = 0;
    cod_info->table_select[1] = 0;
    cod_info->table_select[2] = 0;
    
    {
        if (cod_info->address1 > 0)
        {
            cod_info->table_select[0] = new_choose_table(ix, 0, cod_info->address1);
        }
        if (cod_info->address2 > cod_info->address1)
        {
            cod_info->table_select[1] = new_choose_table(ix, cod_info->address1, cod_info->address2);
        }
        if (cod_info->big_values<<1 > cod_info->address2)
        {
            cod_info->table_select[2] = new_choose_table(ix, cod_info->address2, cod_info->big_values<<1);
        }
    }
}

void subdivide(gr_info *cod_info) 
{
    int scfb_anz = 0;
    int bigvalues_region;
    
    if (!cod_info->big_values) 
    {
        cod_info->region0_count = 0;
        cod_info->region1_count = 0;
    }
    else 
    {
        bigvalues_region = 2*cod_info->big_values;

        int thiscount, index;
        while (sfBandIndex[scfb_anz] < bigvalues_region)
        {
            scfb_anz++;
        }
        cod_info->region0_count = subdv_table[scfb_anz].region0_count;
        thiscount = cod_info->region0_count;
        index = thiscount + 1;
        while (thiscount && (sfBandIndex[index] > bigvalues_region)) 
        {
            thiscount--;
            index--;
        }
        cod_info->region0_count = thiscount;

        cod_info->region1_count = subdv_table[scfb_anz].region1_count;
        index = cod_info->region0_count + cod_info->region1_count + 2;
        thiscount = cod_info->region1_count;
        while (thiscount && (sfBandIndex[index] > bigvalues_region)) 
        {
            thiscount--;
            index--;
        }
        cod_info->region1_count = thiscount;
        cod_info->address1 = sfBandIndex[cod_info->region0_count+1];
        cod_info->address2 = sfBandIndex[cod_info->region0_count
                                                    + cod_info->region1_count+2];
        cod_info->address3 = bigvalues_region;
    }
}

int count1_bitcount(int *ix, gr_info *cod_info) 
{
    int p, i, k;
    int v, w, x, y, signbits;
    int sum0 = 0,
        sum1 = 0;

    for (i=cod_info->big_values<<1, k=0; k<cod_info->count1; i+=4, k++) 
    {   
        v = abs(ix[i]);
        w = abs(ix[i+1]);
        x = abs(ix[i+2]);
        y = abs(ix[i+3]);

        p = v + (w<<1) + (x<<2) + (y<<3);
        
        signbits = 0;
        if (v!=0) 
        {
            signbits++;
        }
        if (w!=0) 
        {
            signbits++;
        }        
        if (x!=0) 
        {
            signbits++;
        }        
        if (y!=0) 
        {
            signbits++;
        }
        sum0 += signbits;
        sum1 += signbits;

        sum0 += ht[32].hlen[p];
        sum1 += ht[33].hlen[p];
    }

    if (sum0<sum1) 
    {
        cod_info->count1table_select = 0;
        return sum0;
    }
    else 
    {
        cod_info->count1table_select = 1;
        return sum1;
    }
}

void calc_runlen(int ix[576], gr_info *cod_info) 
{
    int i;
    for (i = 576; i > 1; i -= 2) 
    {
        if (ix[i-1] || ix[i-2]) 
        {
            break;
        }
    }
    cod_info->count1 = 0 ;
    for (; i > 3; i -= 4) 
    {
        if (abs(ix[i-1]) <= 1 && abs(ix[i-2]) <= 1 && abs(ix[i-3]) <= 1 && abs(ix[i-4]) <= 1)
        {
            cod_info->count1++;
        }
        else
        { 
            break;
        }
    }
    cod_info->big_values = i>>1;
}

void quantize(double xrs[576], int ix[576], gr_info *cod_info ) 
{
    double step;
    double dbl;
    double q;
    step = pow(2.0, (cod_info->quantizerStepSize)/4.0);
    for(int i=0; i<576; i++) 
    {
        dbl = fabs(xrs[i])/step;
        if (dbl>1.0 && dbl<1000000.0) 
        {
            q= pow(dbl, 0.75);
            ix[i]=(int)q;
        }
        else 
        {
        ix[i]=0; //1
        }
    }
}

int count_bits(int ix[576], gr_info *cod_info) 
{
    int bits;

    calc_runlen(ix, cod_info);
     
    if (ix_max(ix, 0, 576) > 8206)
    {
        return 100000;         
    }
    bits = count1_bitcount(ix, cod_info); 
    subdivide(cod_info);
    bigv_tab_select(ix, cod_info); 
    bits += bigv_bitcount(ix, cod_info);
    return bits;
}

int bin_search_StepSize(int desired_rate, double start,
                        double *xrs, int *ix, gr_info *cod_info) 
{
    int top, bot, next, last;
    int bit;

    top  = start;
    next = start;
    bot  = 200;

    do 
    {
        last = next;
        next = ((long)(top+bot) >> 1);
        cod_info->quantizerStepSize = next;
        quantize(xrs, ix, cod_info);
        bit = count_bits(ix, cod_info);
        if (bit > desired_rate)
        {
            top = next;
        }
        else
        {
            bot = next;
        }
    } 
    while((bit!=desired_rate) && abs(last-next)>1);
    return next;
}

int inner_loop(double xrs[576], int ix[576], gr_info *cod_info) 
{
    int bits;
    cod_info->quantizerStepSize += 1.0;
    quantize(xrs, ix, cod_info);
    calc_runlen(ix, cod_info);                  
    bits = count1_bitcount(ix, cod_info);
    subdivide(cod_info);
    bigv_tab_select(ix, cod_info);                
    bits += bigv_bitcount(ix, cod_info);

    while((bits>1588) || (ix_max(ix, 0, 576) > 8206))
    {
        cod_info->quantizerStepSize += 1.0;
        quantize(xrs, ix, cod_info);
        calc_runlen(ix, cod_info);                  
        bits = count1_bitcount(ix, cod_info);
        subdivide(cod_info);
        bigv_tab_select(ix, cod_info);                
        bits += bigv_bitcount(ix, cod_info);
    } 
    return bits;
}

int outer_loop(gr) 
{
    unsigned part23len;
    gr_info *cod_info = &side_info.gr[gr];
    double *xrs = (double *) &(mdct_freq[gr][0]); 
    int *ix = (int *) &(enc[gr][0]);
    bin_search_StepSize(1588, cod_info->quantizerStepSize, xrs, ix, cod_info); 
    part23len = inner_loop(xrs, ix, cod_info);
    return part23len;
}

int quantanf_init(double xr[576]) 
{
    int tp = 0;
    double system_const = 8;
    double sfm = 0, sum1 = 0, sum2 = 0;
    for (int i=0; i<576; i++) 
    {
        if (xr[i]) 
        {   
            double tpd = xr[i]*xr[i];
            sum1 += log(tpd);
            sum2 += tpd;
        }
    }
    if (sum2) 
    {
        sfm = exp(sum1/576)/(sum2/576);
        tp  = (int)(system_const*log(sfm));
        if (tp<-100)
        {
            tp = -100;
        }
    }
    return (tp-70);
}

void code_side_info() 
{
    gr_info* cod_info;     
    unsigned char side_i[18];
    bit_pos(0);
    bit_set(9, side_info.main_data_begin, side_i); 
    bit_set(5, side_info.private_bits, side_i);
    for (int i=0; i<4; i++)
    {
        bit_set(1, (unsigned)side_info.scfsi[i], side_i);
    }
    for (int gr=0; gr<2; gr++) 
    {
        cod_info = (gr_info*) &(side_info.gr[gr]);     
        bit_set(12, cod_info->part2_3_length, side_i);
        bit_set(9, cod_info->big_values, side_i);
        bit_set(8, cod_info->global_gain, side_i);
        bit_set(4, cod_info->scalefac_compress, side_i);
        bit_set(1, cod_info->window_switching_flag, side_i);
        for (int i=0; i<3; i++)
            bit_set(5, cod_info->table_select[i], side_i);
        bit_set(4, cod_info->region0_count, side_i);
        bit_set(3, cod_info->region1_count, side_i);
        bit_set(1, cod_info->preflag, side_i);
        bit_set(1, cod_info->scalefac_scale, side_i);
        bit_set(1, cod_info->count1table_select, side_i);
    } 
    fwrite(side_i, 1, 17, mp3_file); 
}

void iteration_loop()
{
    gr_info *cod_info;

    side_info.main_data_begin = 0;
    side_info.private_bits = 0;

    for(int gr=0; gr<2; gr++) 
    {
        cod_info = (gr_info*) &(side_info.gr[gr]);     
        for (int i=0; i<4; i++) 
        { 
            cod_info->slen[i] = 0;
        }
        cod_info->part2_3_length    = 0;
        cod_info->big_values        = 0;
        cod_info->count1            = 0;
        cod_info->scalefac_compress = 0;
        cod_info->window_switching_flag=0;
        cod_info->table_select[0]   = 0;
        cod_info->table_select[1]   = 0;
        cod_info->table_select[2]   = 0;
        cod_info->region0_count     = 0;
        cod_info->region1_count     = 0;
        cod_info->part2_length      = 0;
        cod_info->preflag           = 0;
        cod_info->scalefac_scale    = 0;
        cod_info->quantizerStepSize = 0.0;
        cod_info->count1table_select= 0;
        cod_info->quantizerStepSize = (double) quantanf_init(mdct_freq[gr]);
        cod_info->part2_3_length  = outer_loop(gr);
        cod_info->global_gain = cod_info->quantizerStepSize+210;
        for (int i=0; i<576; i++) 
        {
            if (mdct_freq[gr][i]<0)
            { 
                enc[gr][i] *=-1;
            }
        }
    }
    code_side_info();
}

static void mdct(double in[36], double out[18]) 
{
    for(int m=0; m<18; m++ ) 
    {
        out[m]= win[35] * in[35] * cos_l[m][35];
        for(int k=0; k<35; k++)
        {
            out[m] += win[k] * in[k] * cos_l[m][k];
        }
    }
}

void mdct_sub() 
{
    double mdct_enc[2][32][18];
    double mdct_in[36];
    
    for(int gr=0; gr<2; gr++) 
    {
        for (int k=0; k<32; k++)
        {
            for (int j=0; j<18; j++)
            {
                mdct_enc[gr][k][j] = mdct_freq[gr][k*18+j];
            }
        }
        for(int band=0; band<32; band++) 
        {
            for(int k=0; k<18; k++)
            {
                if((band&1) && (k&1))
                {
                    sb_sample[gr+1][k][band] *= -1.0;
                }
            }
        }
        for(int band=0; band<32; band++) 
        {
            for(int k=0; k<18; k++) 
            {
                mdct_in[k] = sb_sample[gr][k][band];
                mdct_in[k+18] = sb_sample[gr+1][k][band];
            }
            mdct(mdct_in, &mdct_enc[gr][band][0]);
        }
        for(int band=0; band<31; band++)
        {
            for(int k=0; k<8; k++ ) 
            {
                double bu, bd;
                bu = mdct_enc[gr][band][17-k] * cs[k] + mdct_enc[gr][band+1][k] * ca[k];
                bd = mdct_enc[gr][band+1][k] * cs[k] - mdct_enc[gr][band][17-k] * ca[k];
                mdct_enc[gr][band][17-k] = bu;
                mdct_enc[gr][band+1][k]  = bd;
            }
        }
        for (int k=0; k<32; k++)
        {
            for (int j=0; j<18; j++)
            {
                mdct_freq[gr][k*18+j] = mdct_enc[gr][k][j];
            }
        }
    }
    for(int j=0; j<18; j++) 
    {
        for(int band=0; band<32; band++)
        {
            sb_sample[0][j][band] = sb_sample[2][j][band];
        }
    }
}

void filter_subband(double sb_s[32]) 
{
    double y[64];
    for (int i=0; i<64; i++ ) 
    {
        y[i] = 0;
        for (int j=0; j<8; j++) 
        {
            y[i] += win_que[i+(j<<6)];
        }
    }
    for (int i=0; i<32; i++) 
    {
        sb_s[i] = 0;
        for (int j=0; j<64; j++) 
        {
            sb_s[i] += filter[i][j] * y[j];
        }
    }
}

void window_subband(short buff[32]) 
{
    for (int i=0; i<480; i++)
    { 
        now_que[511-i] = now_que[479-i];
    }
    for (int i=0; i<32;i++)
    {
        now_que[i] = (double) buff[31-i]/32768.0;
    }
    for (int i=0; i<512; i++) 
    {
        win_que[i] = now_que[i]*enwindow[i];
    }
}

void polyph() 
{
    for (int gr=0; gr<2; gr++) 
    {
        for(int i=0; i<18; i++) 
        {
            window_subband(&buffer[gr*576+i*32]);
            filter_subband(sb_sample[gr+1][i]);
        }
    }
}

int wave_get()
{
    int samples_read;

    samples_read = fread(buffer, sizeof(short), 1152, wav_file);
    
    if (samples_read == 1152)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void compress() 
{
    unsigned char head[4] = {0xff, 0xfb, 0x92, 0xc4};
    
    subband_initialise();
    mdct_initialise();
    
    while (wave_get()) 
    {
        polyph();
        mdct_sub();
        fwrite(head, 1, 4, mp3_file); 
        iteration_loop();
        huff();
    }    
}

int main(int argc, char** argv) 
{
    wav_file = fopen(argv[1], "rb");
    fseek(wav_file, 44, SEEK_SET);
    mp3_file = fopen(argv[2], "wb");
    compress();
    fclose(mp3_file);
    fclose(wav_file);
    return 0;
} 

