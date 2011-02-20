#ifndef SSE_H_INCLUDED
#define SSE_H_INCLUDED

#include <xmmintrin.h>
#include <smmintrin.h>

#define loadps(mem)					_mm_load_ps((const float * const)(mem))
#define loadss(mem)					_mm_load_ss((const float * const)(mem))
#define storess(ss,mem)				_mm_store_ss((float * const)(mem),(ss))
#define storeps(ps,mem)				_mm_store_ps((float * const)(mem),(ps))
#define minss						_mm_min_ss
#define maxss						_mm_max_ss
#define minps						_mm_min_ps
#define maxps						_mm_max_ps
#define mulps						_mm_mul_ps
#define mulss						_mm_mul_ss
#define subps						_mm_sub_ps
#define subss						_mm_sub_ss
#define shuffleps					_mm_shuffle_ps
#define rotatelps(ps)				shuffleps((ps),(ps), 0x39)		// a,b,c,d -> b,c,d,a
#define muxhps(low,high)			_mm_movehl_ps((low),(high))		// low{a,b,c,d}|high{e,f,g,h} = {c,d,g,h}
#define dotps						_mm_dp_ps
#define addps						_mm_add_ps
#define addss						_mm_add_ss
#define xorps						_mm_xor_ps
#define movemaskps					_mm_movemask_ps
#define replicate(val)				{val, val, val, val}
#define cmplessps					_mm_cmplt_ps
#define cmpgreatps					_mm_cmpgt_ps
#define cmpleqps					_mm_cmple_ps
#define cmpgteqps					_mm_cmpge_ps
#define cmplessss_i					_mm_comilt_ss
#define setSSE(val)					_mm_set1_ps(val)
#define setZero						_mm_setzero_ps()
#define blendps						_mm_blend_ps


static const float flt_plus_inf = -logf(0);							// let's keep C and C++ compilers happy.
static const float _MM_ALIGN16
	ps_cst_plus_inf[4]	= {  flt_plus_inf,  flt_plus_inf,  flt_plus_inf,  flt_plus_inf },
	ps_cst_minus_inf[4]	= { -flt_plus_inf, -flt_plus_inf, -flt_plus_inf, -flt_plus_inf },
	SSE_invertVec[4]	= {-1,-1,-1,1},
	SSE_invertVec0[4]	= {-1,-1,-1,0};

static const __m128 _zerosps = setZero;
static const __m128 _onesps = setSSE(1);
static const __m128 _negonesps = setSSE(-1);

static __forceinline __m128 hMax(const __m128 & val)
{
	__m128 xmm0(val);
	__m128 xmm1(muxhps( xmm0, xmm0 ));
	xmm0 = maxps( xmm0, xmm1 );
	xmm1 = shuffleps( xmm0, xmm0, _MM_SHUFFLE(1,1,1,1) );
	return maxss( xmm0, xmm1 );
}

static __forceinline __m128 hMin(const __m128 & val)
{
	__m128 xmm0(val);
	__m128 xmm1(shuffleps( xmm0, xmm0, _MM_SHUFFLE(2,2,2,2)));
	xmm0 = minps( xmm0, xmm1 );
	xmm1 = shuffleps( xmm0, xmm0, _MM_SHUFFLE(1,1,1,1) );
	return minss( xmm0, xmm1 );
}

static __forceinline __m128 recipss(const __m128 & val)
// Reciprocal of the lowest float in the SSE register.
{
	__m128 xmm0(_mm_rcp_ss(val));
	return subss(mulss(setSSE(2.0f), xmm0), mulss(val, mulss(xmm0, xmm0))); // One step of Newton-Raphson, otherwise accuracy is too low.
}

static __forceinline void recipss(const __m128 & val, float & out)
// Reciprocal of the lowest float in the SSE register.
{
	__m128 xmm0(_mm_rcp_ss(val));
	storess(subss(mulss(setSSE(2.0f), xmm0), mulss(val, mulss(xmm0, xmm0))), &out); // One step of Newton-Raphson, otherwise accuracy is too low.
}

static __forceinline __m128 recipps(const __m128 & val)
// Reciprocal of all floats in the SSE register.
{
	__m128 xmm0(_mm_rcp_ps(val));
	return subps(mulps(setSSE(2.0f), xmm0), mulps(val, mulps(xmm0, xmm0))); // One step of Newton-Raphson, otherwise accuracy is too low.
}

static __forceinline __m128 fastrsqrtps(const __m128 & val)
// Fast reciprocal square root of packed scalars. Also uses Newton-Raphson to increase accuracy
{
	const __m128 approx = _mm_rsqrt_ps(val);
	const __m128 muls = mulps(mulps(val, approx), approx);
	return mulps(mulps(setSSE(0.5f), approx), subps(setSSE(3), muls) );
}

static __forceinline void fastrsqrtss(const __m128 & val, float & out)
{
	const __m128 approx = _mm_rsqrt_ss(val);
	const __m128 muls = mulss(mulss(val, approx), approx);
	storess(mulss(mulss(setSSE(0.5f), approx), subss(setSSE(3), muls) ), &out);
}

static __forceinline void SoACross(const __m128 & X1, const __m128 & Y1, const __m128 & Z1, const __m128 & X2, const __m128 & Y2, const __m128 & Z2, __m128 & XR, __m128 & YR, __m128 & ZR)
// Structure-of-Array
{
	XR = subps(mulps(Y1,Z2), mulps(Z1,Y2));
	YR = mulps(setSSE(-1), subps(mulps(X1,Z2), mulps(Z1,X2)));
	ZR = subps(mulps(X1,Y2), mulps(Y1,X2));
}

static __forceinline __m128 SoADot(const __m128 & X1, const __m128 & Y1, const __m128 & Z1, const __m128 & X2, const __m128 & Y2, const __m128 & Z2)
{
	return addps(mulps(X1, X2), addps(mulps(Y1, Y2), mulps(Z1, Z2)));
}

#endif