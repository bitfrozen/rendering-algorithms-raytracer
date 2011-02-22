#ifndef CSE168_BLINN_H_INCLUDED
#define CSE168_BLINN_H_INCLUDED

#include "Material.h"
#include "mtrand.h"
#include "Miro.h"
#include "SSE.h"

class Blinn : public Material
{
public:
    Blinn(const Vector3 & kd = Vector3(1.f),
			const Vector3 & ka = Vector3(0.f),
			const Vector3 & ks = Vector3(1.f),
			const Vector3 & kt = Vector3(0.f),
			float ior = 1.5f,
			float specExp = 1.0,
			float specAmt = 1.0,
			float reflectAmt = 0.0,
			float refractAmt = 0.0);
    virtual ~Blinn();

    const Vector3 & kd() const {return m_kd;}
    const Vector3 & ka() const {return m_ka;}
	const Vector3 & ks() const {return m_ks;}
	const Vector3 & kt() const {return m_kt;}
	const float ior() const {return m_ior;}
	const float specExp() const {return m_specExp;}
	const float specAmt() const {return m_specAmt;}
	const float reflectAmt() const {return m_reflectAmt;}
	const float refractAmt() const {return m_refractAmt;}

    void setKd(const Vector3 & kd) {m_kd = kd;}
    void setKa(const Vector3 & ka) {m_ka = ka;}
	void setKs(const Vector3 & ks) {m_ks = ks;}
	void setKt(const Vector3 & kt) {m_kt = kt;}
	void setIor(const float ior) {m_ior = ior;}
	void setSpecExp(const float specExp) {m_specExp = specExp;}
	void setSpecAmt(const float specAmt) {m_specAmt = specAmt;}
	void setReflectAmt(const float reflectAmt) {m_reflectAmt = reflectAmt;}
	void setRefractAmt(const float refractAmt) {m_refractAmt = refractAmt;}
	
	void setLightEmittedIntensity(float le) { m_lightEmitted = le; }
	void setLightEmittedColor(const Vector3 & le) { m_Le = le; }
    virtual void preCalc() {}
    
    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;
	static int randsIdx;
	static float rands[1000000];
	static void genRands();

protected:
    Vector3 m_kd;			// Diffuse Color
    Vector3 m_ka;			// Ambient Color
	Vector3 m_ks;			// Specular / Reflection Color
	Vector3 m_kt;			// Transmittance (Refraction) Color
	float m_ior;			// Index of Refraction
	float m_specExp;		// Specular exponent
	float m_specAmt;		// Specular component weight (only for "hilights")
	float m_reflectAmt;		// Reflection amount. 1.0 -> 100% reflectve. Weighted using fresnel approximation.
	float m_refractAmt;		// Refraction amount. This weights the refraction amount prescribed by the fresnel approximation.
	float m_lightEmitted;	//power of light emitted
	Vector3 m_Le;			//color of the emitted light
};

#endif // CSE168_LAMBERT_H_INCLUDED