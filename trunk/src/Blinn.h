#ifndef CSE168_BLINN_H_INCLUDED
#define CSE168_BLINN_H_INCLUDED

#include "Material.h"

class Blinn : public Material
{
public:
    Blinn(const Vector3 & kd = Vector3(1),
			const Vector3 & ka = Vector3(0),
			const Vector3 & ks = Vector3(1),
			const Vector3 & kt = Vector3(0),
			float ior = 1.5,
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

    virtual void preCalc() {}
    
    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;
protected:
    Vector3 m_kd;
    Vector3 m_ka;
	Vector3 m_ks;
	Vector3 m_kt;
	float m_ior;
	float m_specExp;
	float m_specAmt;
	float m_reflectAmt;
	float m_refractAmt;
};

#endif // CSE168_LAMBERT_H_INCLUDED
