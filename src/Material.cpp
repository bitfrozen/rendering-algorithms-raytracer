#include "Material.h"

Material::Material() : m_envExposure(1.0f), m_envMap(NULL), m_texture(NULL)
{
}

Material::~Material()
{
}

Vector3
Material::shade(const Ray&, const HitInfo&, const Scene&) const
{
    return Vector3(1.0f, 1.0f, 1.0f);
}
