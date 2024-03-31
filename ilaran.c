// WELL Algorithm for generating random integers
// It performs much better than rand() but is still simple
// by Yiping Cheng, Beijing Jiaotong University, China
// Email: ypcheng@bjtu.edu.cn

typedef unsigned int UINT;

// initialize state to random bits
static UINT state[16] =
{
    0xED441570, 0x2BBB76C8, 0x73276792, 0x66739C9C,
    0x6189D44A, 0x2C1BAC00, 0x9DEC31E8, 0xB0100E59,
    0xCDDF7074, 0x42C707D0, 0x6AAE869C, 0x2F159D0D,
    0x24AA6EEB, 0x3786C1C7, 0xB705D8A4, 0xC89578C0
};

// init should also reset this to 0
static int index = 0;

void ILASEED(UINT sid)
{
    UINT* ptr = state;
    int i;

    for (i = 1; i < 9; i++) {
        ptr[0] = sid * i;
        ptr[1] = ~ptr[0];

        ptr += 2;
    }

    index = 0;
}


// generates 32-bit random integer
UINT ILARAN()
{
    UINT a, b, c, d;

    a = state[index];
    c = state[(index + 13) & 15];
    b = a ^ c ^ (a << 16) ^ (c << 15);
    c = state[(index + 9) & 15];
    c ^= (c >> 11);
    a = state[index] = b ^ c;
    d = a ^ ((a << 5) & 0xDA442D24U);
    index = (index + 15) & 15;
    a = state[index];

    return state[index] = a ^ b ^ d ^ (a << 2) ^ (b << 18) ^ (c << 28);
}
