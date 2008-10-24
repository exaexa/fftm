
/*
 * FFT long integer operations
 * by Mirek Kratochvil ( [exa] ) 2008
 *
 * You may only read, use, modify and redistribute this file,
 * this all under terms of GNU GPLv3.
 *
 * Please see www.gnu.org for details
 *
 *
 * How to compile
 *
 * g++ -o fftc fftc.cpp
 *
 * you might also need to pass -lm argument, if your system doesn't
 * link it by default.
 */


#include <stdio.h>
#include <string.h>
#include <math.h>

class complex
/*
 * small class for complex convenience.
 * just a small percent of implemented operations is used here, but I'm just
 * somewhat lazy to clean it out.
 */
{

public:

	union{
		struct{
			double x, y;
		};
		double v[2];
	};

	complex (double a = 0, double b = 0) : x (a), y (b)
	{}

	complex (const complex & a)
	{
		x = a.x;
		y = a.y;
	}

	complex operator+ (const complex & a) const
	{
		return complex (a.x + x, a.y + y);
	}

	complex operator+()
	{
		return *this;
	}

	complex operator- (const complex & a) const
	{
		return complex (x -a.x, y - a.y);
	}

	complex operator-()
	{
		return complex (-x, -y);
	}

	complex operator* (double a) const
	{
		return complex (a*x, a*y);
	}

	friend complex operator* (double a, const complex & c)
	{
		return complex (a* (c.x), a* (c.y) );
	}

	complex operator/ (double a) const
	{
		return complex (x / a, y / a);
	}

	complex operator* (const complex& a) const
	{
		return complex ( (x*a.x) - (y*a.y) , (x*a.y) + (y*a.x) );
	}

	complex operator/ (const complex& a) const
	{
		return (*this*complex (a.x, -a.y) ) / (a.x*a.x + a.y*a.y);
	}

	const complex& operator+= (const complex& a)
	{
		x += a.x;
		y += a.y;
		return *this;
	}

	const complex& operator-= (const complex& a)
	{
		x -= a.x;
		y -= a.y;
		return *this;
	}

	const complex operator*= (const double a)
	{
		(*this) = (*this) * a;
		return *this;
	}
	const complex operator*= (const complex& a)
	{
		(*this) = (*this) * a;
		return *this;
	}
	const complex operator/= (const double a)
	{
		(*this) = (*this) / a;
		return *this;
	}
	const complex operator/= (const complex& a)
	{
		(*this) = (*this) / a;
		return *this;
	}

	double length() const
	{
		return (double) pow ( (x*x) + (y*y), 0.5f);
	}
	complex unit() const
	{
		return (*this) / length();
	}
	const complex& normalize()
	{
		return *this = unit();
	}

	complex operator| (const double& size) const
	{
		return size*unit();
	}
	const complex & operator|= (const double& size)
	{
		(*this) = size * unit();
		return *this;
	}

	double operator% (const complex& a)
	{  //dot product
		return x*a.x + y*a.y;
	}

	bool operator== (const complex& a) const
	{
		return (x == a.x) && (y == a.y);
	}
	bool operator!= (const complex& a) const
	{
		return (x != a.x) || (y != a.y);
	}
};


/*
 * Compute FFT in direction forward (polynome->values) or reverse (<-).
 * x is an array of complex parameters exactly 2^m long.
 *
 * eg. Thanks to Paul Bourke, fourier code is greatly inspired by his
 * optimalizations. Especially the integer magic. kthx.
 */

void FFT (bool forward, int m, complex x[])
{
	int i, i1, i2, j, k, l, l1, l2, n;
	complex t, u, c;

	//compute real size
	n = 1;
	for (i = 0; i < m; i++)
		n <<= 1;

	//swap polynomials so they are matched in pairs
	i2 = n >> 1;
	j = 0;

	for (i = 0; i < n - 1 ; i++) {
		if (i < j) {
			t = x[i];
			x[i] = x[j];
			x[j] = t;
		}

		k = i2;

		while (k <= j) {
			j -= k;
			k >>= 1;
		}

		j += k;
	}

	//FFT
	c = complex (-1, 0);
	l2 = 1;
	for (l = 0; l < m; l++) {
		l1 = l2;
		l2 <<= 1;

		u = complex (1, 0);

		for (j = 0; j < l1; j++) {
			for (i = j; i < n; i += l2) {
				i1 = i + l1;
				t = u * x[i1];
				x[i1] = x[i] - t;
				x[i] += t;
			}

			u = u * c;
		}

		c = complex (sqrt ( (1 + c.x) / 2), sqrt ( (1 - c.x) / 2) );
		if (!forward) c.y *= -1;
	}

	//divide by argument number. backward fft needs this.
	if (!forward) {
		for (i = 0; i < n; i++)
			x[i] /= n;
	}
}

/*
 * Here we come with a plenty of help functions. Basically we load numbers
 * into a linked list, then measure the size of it, create a filed of complex
 * numbers of next_power_of_2 size, FFT it, multiply it, FFT it back, truncate
 * the noise that might be generated, then print it out.
 */


//conversion to the powers of 2
int next_p2 (int a)
{
	int r = 1;
	while ( (r < a) && r) r <<= 1;
	return r;
}

int which_p2 (int a)
{
	int r=0;
	while(a!=1) {a>>=1; ++r;}
	return r;
}

//linked lists
typedef struct number_list {
	int a;
	struct number_list *next;
} *nlist;

void printlist(nlist a)
{
	while(a) {printf("%d -> ",a->a); a=a->next;}
	printf("null\n");
}

void newlist(nlist*l)
{
	*l=0;
}

void push(nlist *l,int val)
{
	nlist nl=new struct number_list;
	nl->a=val;
	nl->next=*l;
	*l=nl;
}

bool pop(nlist *l,int*val)
{
	if(!(*l))return false;
	if(val) *val=(*l)->a;
	nlist t=*l;
	*l=t->next;
	delete t;
}

void clean(nlist*l)
{
	while(pop(l,0));
}

void revlist(nlist*l)
{
	nlist last=0,t;
	while(*l){
		t=(*l);
		(*l)=t->next;
		t->next=last;
		last=t;
	}
	*l=last;
}

//read arbitrarily long number and save it to linked list
int read_number_to_list(FILE*f, nlist*l)
{
	char c;
	int n=0;
	newlist(l);
	while(((c=fgetc(f))!='\n')&&(!feof(f)))
		if(((c>='0')&&(c<='9'))
			&&((c!='0')||(n>0))) //starting zeroes (00023 => 23)
		{
			push(l,c-'0');
			++n;
		}
	return n;
}

//print the number in list
void print_list_as_number(FILE*f, nlist l)
{
	if(!l)fprintf(f,"%d",0);
	revlist(&l);
	
	for(;l;l=l->next)fprintf(f,"%d",l->a);

	revlist(&l);
	fprintf(f,"\n");
}

//conversion list->field
void list_to_fld(complex*c, nlist a, int len, int maxlen)
{
	int i;
	for(i=0;a&&(i<len);++i,(a=a->next)) c[i]=complex(a->a,0);
	for(;i<maxlen;++i) c[i]=complex(0,0);
}


//Note that field->list conversion has to be messed up a little with special
//rounding rules/tweaks.

int do_round(float f, float percision)
{
	if(f<0)return 0;
	if(f<percision)return 0;
	if(f<1)return 1;
	return (int) floorf (f);
}

void fld_round_and_move_overflows(complex*c, int maxlen, float percision)
{
	int i,t=0;
	for(i=0;i<maxlen;++i){
		t+=do_round(c[i].x,percision);
		c[i].x=t%10;
		t/=10;
	}
}

int fld_to_list(complex*c, nlist*a, int maxlen, float percision)
{
	int i,t;
	fld_round_and_move_overflows(c,maxlen,percision);
	for(i=maxlen-1;(i>=0)&&(!(int)(c[i].x));--i); //strip zeroes
	for(;i>=0;--i)
		push(a,(int)(c[i].x));
}

//the multiply operation itself.
int do_mul(nlist a1, nlist a2, int len1, int len2, nlist* result)
{
	newlist(result);
	if(!(len1&&len2)) return 0; //GAH ZERO MULTIPLY!
	complex *n1,*n2;
	int len=next_p2(len1+len2);
	int lenp2=which_p2(len);

	n1=new complex[len];
	n2=new complex[len];
	
	list_to_fld(n1,a1,len1,len);
	list_to_fld(n2,a2,len2,len);
	
	FFT(true, lenp2, n1);
	FFT(true, lenp2, n2);

	for(int i=0;i<len;++i) n1[i]*=n2[i];

	FFT(false, lenp2, n1);

	int r=fld_to_list(n1,result,len,0.05);

	delete [] n1; delete [] n2;
	return r;
}

char helpstr[]=
"\
Compute an operation on two numbers using FFT\n\
\n\
    USAGE:\n\
\n\
Program reads two decimal integral numbers from the specified file,\n\
Result of multiplying will be output to another file.\n\
Default files are stdin and stdout.\n\
\n\
    Options:\n\
\n\
-i <file>  input file\n\
-o <file>  output file (will be overwritten)\n\
\n";



int main(int argc, char**argv)
{
	nlist a1,a2,r;
	int l1,l2,lr;

	FILE*i=0;
	FILE*o=0;

	int mode=0;
	int printinfo=0;


	for(int arg=1;arg<argc;++arg) //argument parsing
		switch(mode){
		case 1: //in filename
			if(i)fclose(i);
			i=fopen(argv[arg],"r");
			mode=0;
			break;
		case 2: //out filename
			if(o)fclose(o);
			o=fopen(argv[arg],"w");
			mode=0;
			break;
		default:
			if(!strcmp(argv[arg],"-i"))
				mode=1;
			else if(!strcmp(argv[arg],"-o"))
				mode=2;
			else if(!strcmp(argv[arg],"--help")) printinfo|=1;
			else if(!strcmp(argv[arg],"--version")) printinfo|=2;
			else printf("invalid argument `%s' not parsed!\n",
				argv[arg]);
		}
		
	if(printinfo&2)printf("FFT Calc by [exa] version 1.0.0.0 pwn!\n");
	if(printinfo&1)printf(helpstr);
	if(printinfo) goto skip;

	l1=read_number_to_list(i?i:stdin,&a1);
	l2=read_number_to_list(i?i:stdin,&a2);
	lr=do_mul(a1,a2,l1,l2,&r);

	print_list_as_number(o?o:stdout,r);

	clean(&a1);
	clean(&a2);
	clean(&r);

skip:	
	if(i)fclose(i);
	if(o)fclose(o);

	return 0;
}

