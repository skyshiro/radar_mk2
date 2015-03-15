#include "msp430.h"
#include <stdlib.h>
#include <stdint.h>

//Latch port allows user to change the pin used to latch CS
void DAC_cipher(int amplitude, int latch_port);
void DAC_setup();

void SetClock(unsigned int clkSpeed);
void SetVcoreUp (unsigned int level);


#define DAC_SAMPLE 250
#define DAC_CONST 4095/DAC_SAMPLE

/* Select the global Q value and include the Qmath header file. */
#define GLOBAL_Q    12
#include "QmathLib.h"

/* Specify the sample size and sample frequency. */
#define FFT_SAMPLE              256      // <= 256, power of 2
#define SAMPLE_FREQUENCY        25000    // <= 16384
#define PI						3.14259

/* Access the real and imaginary parts of an index into a complex array. */
#define RE(x)           (((x)<<1)+0)    // access real part of index
#define IM(x)           (((x)<<1)+1)    // access imaginary part of index

/*
 * Input and result buffers. These can be viewed in memory or printed by
 * defining ALLOW_PRINTF.
 */
_q qInput[FFT_SAMPLE*2];                   // Input buffer of complex values
_q qMag[FFT_SAMPLE/2];                     // Magnitude of each frequency result
_q qPhase[FFT_SAMPLE/2];                   // Phase of each frequency result

extern void cFFT(_q *input, int16_t n);

int adcFlag = 0;
int DAC_flag;

int sample_count;
int offset = 1150;  //934
int amplitude = 1150;
int m;
int big_avg = 0;
int small_avg = 0;


int main(void)
{
    int16_t i, j;                       // loop counters
    int samples[FFT_SAMPLE],k=0;
    int period = 10000; // in us

    /* Disable WDT. */
    WDTCTL = WDTPW + WDTHOLD;
    
    SetClock(16000);

    ADC12CTL0 = ADC12ON;         // Sampling time, ADC12 on ADC12SHT01 +
    ADC12CTL1 = ADC12SHP + ADC12SSEL1;        // Use sampling timer with MCLK
    ADC12CTL2 &= ~ADC12RES1;					// 8bit resolution
    ADC12IE = 0x01;                           // Enable interrupt
    ADC12CTL0 |= ADC12ENC;
    P6SEL |= 0x01;                            // P6.0 ADC option select
    P1DIR |= 0x01;                            // P1.0 output
    __bis_SR_register(GIE);     // LPM0, ADC12_ISR will force exit

    DAC_setup();

	TA0CCTL0 = CCIE;                            // CCR0 interrupt enabled
	TA0CCR0 = period*2/DAC_SAMPLE;
	TA0CTL = TASSEL_2 + MC_1 + TAIE + ID_3;     // SMCLK, up mode, /8 for 2 MHz clock. MC_1 turns on timer

	//Timer resoultion = T_period/(MAX_samples*T_clk)
	//For 1 us clock and T_period in us = period/SAMPLE_MAX

	sample_count = 0;
	DAC_flag = 0;

    while (sample_count < DAC_SAMPLE)
    {
		//DAC_flag is used to determine time to change value of DAC
		if(DAC_flag)
		{
			ADC12CTL0 |= ADC12SC;                   // Start sampling/conversion

			DAC_cipher(DAC_CONST*sample_count, BIT0);

			DAC_flag = 0;
			sample_count++;
		}

		if(adcFlag)
		{
			qInput[sample_count] = ADC12MEM0;
			adcFlag = 0;
		}
    }

    TA0CTL &= ~MC_1;
    TA0CCTL0 &= ~CCIE;

    big_avg = 0;

    for(i=0;i<16;i++)
    {
    	small_avg = 0;

    	for(j=0;j<16;j++)
    	{
    		small_avg += qInput[16*j+i];
    	}

    	big_avg += small_avg/16;
    }

    big_avg /= 16;

    for(i=0;i<FFT_SAMPLE;i++)
    {
    	qInput[i] = qInput[i] - big_avg;
    }
    /*
     * Perform a complex FFT on the input samples. The result is calculated
     * in-place and will be stored in the input buffer.
     */
    cFFT(qInput, FFT_SAMPLE);
    
    /* Calculate the magnitude and phase angle of the results. */
    for (i = 0; i < FFT_SAMPLE/2; i++) {
        qMag[i] = _Qmag(qInput[RE(i)], qInput[IM(i)]);
        qPhase[i] = _Qatan2(qInput[IM(i)], qInput[RE(i)]);
    }
    
    while(1);
    
    return 0;
}

extern void cBitReverse(_q *input, int16_t n);

/*
 * Perform in-place radix-2 DFT of the input signal with size n.
 *
 * This function has been written for any input size up to 256. This function
 * can be optimized by using lookup tables with precomputed twiddle factors for
 * a fixed sized FFT, using Q15 format for the twiddle factors and inlining the
 * multiplication steps with direct access to the MPY32 hardware peripheral.
 */
void cFFT(_q *input, int16_t n)
{
    int16_t s, s_2;                     // step
    uint16_t i, j;                      // loop counters
    _q qTAngle;                         // twiddle factor angle
    _q qTIncrement;                     // twiddle factor increment
    _q qTCos, qTSin;                    // complex components of twiddle factor
    _q qTempR, qTempI;                  // temp result complex pair
    
    /* Bit reverse the order of the inputs. */
    cBitReverse(input, n);
    
    /* Set step to 2 and initialize twiddle angle increment. */
    s = 2;
    s_2 = 1;
    qTIncrement = _Q(-2*PI);
    
    while (s <= n) {
        /* Reset twiddle angle and halve increment factor. */
        qTAngle = 0;
        qTIncrement = _Qdiv2(qTIncrement);
        
        for (i = 0; i < s_2; i++) {
            /* Calculate twiddle factor complex components. */
            qTCos = _Qcos(qTAngle);
            qTSin = _Qsin(qTAngle);
            qTAngle += qTIncrement;
            
            for (j = i; j < n; j += s) {
                /* Multiply complex pairs and scale each stage. */
                qTempR = _Qmpy(qTCos, input[RE(j+s_2)]) - _Qmpy(qTSin, input[IM(j+s_2)]);
                qTempI = _Qmpy(qTSin, input[RE(j+s_2)]) + _Qmpy(qTCos, input[IM(j+s_2)]);
                input[RE(j+s_2)] = _Qdiv2(input[RE(j)] - qTempR);
                input[IM(j+s_2)] = _Qdiv2(input[IM(j)] - qTempI);
                input[RE(j)] = _Qdiv2(input[RE(j)] + qTempR);
                input[IM(j)] = _Qdiv2(input[IM(j)] + qTempI);
            }
        }
        /* Multiply step by 2. */
        s_2 = s;
        s = _Qmpy2(s);
    }
}

/*
 * Perform an in-place bit reversal of the complex input array with size n.
 * Use a look up table to speed up the process. Valid for size of 256 and
 * smaller.
 */
void cBitReverse(_q *input, int16_t n)
{
    uint16_t i, j;                      // loop counters
    int16_t i16BitRev;                  // index bit reversal
    _q qTemp;
    
    extern const uint8_t ui8BitRevLUT[256];
    
    /* In-place bit-reversal. */
    for (i = 0; i < n; i++) {
        i16BitRev = ui8BitRevLUT[i];
        for (j = n; j < 256; j <<= 1) {
            i16BitRev >>= 1;
        }
        if (i < i16BitRev) {
            /* Swap inputs. */
            qTemp = input[RE(i)];
            input[RE(i)] = input[RE(i16BitRev)];
            input[RE(i16BitRev)] = qTemp;
            qTemp = input[IM(i)];
            input[IM(i)] = input[IM(i16BitRev)];
            input[IM(i16BitRev)] = qTemp;
        }
    }
}

/* 8-bit reversal lookup table. */
const uint8_t ui8BitRevLUT[256] = {
    0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0, 
    0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8, 
    0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4, 
    0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC, 
    0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2, 
    0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
    0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6, 
    0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
    0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
    0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9, 
    0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
    0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
    0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3, 
    0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
    0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7, 
    0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
};


#if defined(__TI_COMPILER_VERSION__) || defined(__IAR_SYSTEMS_ICC__)
#pragma vector = ADC12_VECTOR
__interrupt void ADC12_ISR(void)
#elif defined(__GNUC__)
void __attribute__ ((interrupt(ADC12_VECTOR))) ADC12_ISR (void)
#else
#error Compiler not supported!
#endif
{
  switch(__even_in_range(ADC12IV,34))
  {
  case  0: break;                           // Vector  0:  No interrupt
  case  2: break;                           // Vector  2:  ADC overflow
  case  4: break;                           // Vector  4:  ADC timing overflow
  case  6:                                  // Vector  6:  ADC12IFG0
    adcFlag = 1;
    ADC12IFG &= ~BIT0;
    break;
  case  8: break;                           // Vector  8:  ADC12IFG1
  case 10: break;                           // Vector 10:  ADC12IFG2
  case 12: break;                           // Vector 12:  ADC12IFG3
  case 14: break;                           // Vector 14:  ADC12IFG4
  case 16: break;                           // Vector 16:  ADC12IFG5
  case 18: break;                           // Vector 18:  ADC12IFG6
  case 20: break;                           // Vector 20:  ADC12IFG7
  case 22: break;                           // Vector 22:  ADC12IFG8
  case 24: break;                           // Vector 24:  ADC12IFG9
  case 26: break;                           // Vector 26:  ADC12IFG10
  case 28: break;                           // Vector 28:  ADC12IFG11
  case 30: break;                           // Vector 30:  ADC12IFG12
  case 32: break;                           // Vector 32:  ADC12IFG13
  case 34: break;                           // Vector 34:  ADC12IFG14
  default: break;
  }
}

// Timer A0 interrupt service routine
// DAC_flag is raised in ISR to change DAC value
#pragma vector=TIMER0_A0_VECTOR
__interrupt void Timer_A (void)
{
	//int timerVal = TAR;
	TA0CTL &= ~TAIFG;
	DAC_flag = 1;
}

void DAC_setup()
{
	P4DIR |= BIT0;						// Will use P2.0 to activate /CE on the DAC
	P3SEL	= BIT0 + BIT2;	// + BIT4;	// SDI on P3.0 and SCLK on P3.2

	UCB0CTL0 |= UCCKPL + UCMSB + UCMST + /* UCMODE_2 */ + UCSYNC;
	UCB0CTL1 |= UCSSEL_2;	// UCB0 will use SMCLK as the basis for

	//Divides UCB0 clk, I think @ 4 MHz CPU clk it should be fine
	UCB0BR0 |= 0x10;	// (low divider byte)
	UCB0BR1 |= 0x00;	// (high divider byte)
	UCB0CTL1 &= ~UCSWRST;	// **Initialize USCI state machine**
}

void DAC_cipher(int amplitude, int latch_port)
{
	int DAC_code;

	DAC_code = (0x3000)|(amplitude & 0xFFF); //Gain to 1 and DAC on

	P4OUT &= ~latch_port; //Lower CS pin low
	//send first 8 bits of code
	UCB0TXBUF = (DAC_code >> 8);

	//wait for code to be sent
	while(!(UCTXIFG & UCB0IFG));

	//send last 8 bits of code
	UCB0TXBUF = (char)(DAC_code & 0xFF);

	//wait for code to be sent
	while(!(UCTXIFG & UCB0IFG));

	//wait until latch
	_delay_cycles(100);
	//raise CS pin
	P4OUT |= latch_port;

	return;
}

void SetVcoreUp (unsigned int level)
{
  // Open PMM registers for write
  PMMCTL0_H = PMMPW_H;
  // Set SVS/SVM high side new level
  SVSMHCTL = SVSHE + SVSHRVL0 * level + SVMHE + SVSMHRRL0 * level;
  // Set SVM low side to new level
  SVSMLCTL = SVSLE + SVMLE + SVSMLRRL0 * level;
  // Wait till SVM is settled
  while ((PMMIFG & SVSMLDLYIFG) == 0);
  // Clear already set flags
  PMMIFG &= ~(SVMLVLRIFG + SVMLIFG);
  // Set VCore to new level
  PMMCTL0_L = PMMCOREV0 * level;
  // Wait till new level reached
  if ((PMMIFG & SVMLIFG))
    while ((PMMIFG & SVMLVLRIFG) == 0);
  // Set SVS/SVM low side to new level
  SVSMLCTL = SVSLE + SVSLRVL0 * level + SVMLE + SVSMLRRL0 * level;
  // Lock PMM registers for write access
  PMMCTL0_H = 0x00;
}

void SetClock(unsigned int clkSpeed)
{
  // Increase Vcore setting to level3 to support fsystem=25MHz
  // NOTE: Change core voltage one level at a time..
  SetVcoreUp (0x01);
  SetVcoreUp (0x02);
  SetVcoreUp (0x03);

  UCSCTL3 = SELREF_2;                       // Set DCO FLL reference = REFO
  UCSCTL4 |= SELA_2;                        // Set ACLK = REFO

  __bis_SR_register(SCG0);                  // Disable the FLL control loop
  UCSCTL0 = 0x0000;                         // Set lowest possible DCOx, MODx
  UCSCTL1 = DCORSEL_7;                      // Select DCO range 50MHz operation
  UCSCTL2 = FLLD_0 + (clkSpeed/32.768-1);                   // Set DCO Multiplier for 25MHz
                                            // (N + 1) * FLLRef = Fdco
                                            // (762 + 1) * 32768 = 25MHz
                                            // Set FLL Div = fDCOCLK/2
  __bic_SR_register(SCG0);                  // Enable the FLL control loop

  // Worst-case settling time for the DCO when the DCO range bits have been
  // changed is n x 32 x 32 x f_MCLK / f_FLL_reference. See UCS chapter in 5xx
  // UG for optimization.
  // 32 x 32 x 25 MHz / 32,768 Hz ~ 780k MCLK cycles for DCO to settle
  __delay_cycles(782000);

  // Loop until XT1,XT2 & DCO stabilizes - In this case only DCO has to stabilize
  do
  {
    UCSCTL7 &= ~(XT2OFFG + XT1LFOFFG + DCOFFG);
                                            // Clear XT2,XT1,DCO fault flags
    SFRIFG1 &= ~OFIFG;                      // Clear fault flags
  }while (SFRIFG1&OFIFG);                   // Test oscillator fault flag

 }


