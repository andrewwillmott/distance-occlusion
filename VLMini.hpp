#ifndef VL_MINI_H
#define VL_MINI_H

#include <math.h>
#include <stdint.h>

template<class T> inline T sqr(T x) { return x * x; }

template<class T> inline T vl_min(T a, T b) { return a < b ? a : b; }
template<class T> inline T vl_max(T a, T b) { return a > b ? a : b; }

constexpr float vlf_pi = 3.14159265358979323846f;
constexpr float vlf_halfPi = 0.5 * vlf_pi;
constexpr float vlf_twoPi = 2 * vlf_pi;

#define VI(T) T operator[](int i) const { return (&x)[i]; } T& operator[](int i){ return (&x)[i]; }

struct Vec2f { float x; float y;                    Vec2f() {}; Vec2f(float xi, float yi)                     : x(xi), y(yi)               {}; VI(float) };
struct Vec3f { float x; float y; float z;           Vec3f() {}; Vec3f(float xi, float yi, float zi)           : x(xi), y(yi), z(zi)        {}; VI(float) };
struct Vec4f { float x; float y; float z; float w;  Vec4f() {}; Vec4f(float xi, float yi, float zi, float wi) : x(xi), y(yi), z(zi), w(wi) {}; VI(float) };

struct Vec2d { double x; double y;                      Vec2d() {}; Vec2d(double xi, double yi)                       : x(xi), y(yi)               {}; VI(double) };
struct Vec3d { double x; double y; double z;            Vec3d() {}; Vec3d(double xi, double yi, double zi)            : x(xi), y(yi), z(zi)        {}; VI(double) };
struct Vec4d { double x; double y; double z; double w;  Vec4d() {}; Vec4d(double xi, double yi, double zi, double wi) : x(xi), y(yi), z(zi), w(wi) {}; VI(double) };

struct Mat2f { Vec2f x; Vec2f y;                    Mat2f() {}; Mat2f(Vec2f xi, Vec2f yi)                     : x(xi), y(yi)               {}; VI(Vec2f) };
struct Mat3f { Vec3f x; Vec3f y; Vec3f z;           Mat3f() {}; Mat3f(Vec3f xi, Vec3f yi, Vec3f zi)           : x(xi), y(yi), z(zi)        {}; VI(Vec3f) };
struct Mat4f { Vec4f x; Vec4f y; Vec4f z; Vec4f w;  Mat4f() {}; Mat4f(Vec4f xi, Vec4f yi, Vec4f zi, Vec4f wi) : x(xi), y(yi), z(zi), w(wi) {}; VI(Vec4f) };

#undef VI

inline Vec2f  operator+ (Vec2f  v)          { return { +v.x, +v.y }; }
inline Vec2f  operator- (Vec2f  v)          { return { -v.x, -v.y }; }
inline Vec2f  operator+ (Vec2f  a, Vec2f b) { return { a.x + b.x, a.y + b.y }; }
inline Vec2f  operator- (Vec2f  a, Vec2f b) { return { a.x - b.x, a.y - b.y }; }
inline Vec2f  operator* (Vec2f  a, Vec2f b) { return { a.x * b.x, a.y * b.y }; }
inline Vec2f  operator/ (Vec2f  a, Vec2f b) { return { a.x / b.x, a.y / b.y }; }
inline Vec2f  operator* (float  s, Vec2f a) { return { s * a.x, s * a.y }; }
inline Vec2f  operator* (Vec2f  a, float s) { return { s * a.x, s * a.y }; }
inline Vec2f  operator/ (float  s, Vec2f a) { return { s / a.x, s / a.y }; }
inline Vec2f  operator/ (Vec2f  a, float s) { return { s / a.x, s / a.y }; }
inline Vec2f& operator+=(Vec2f& a, Vec2f b) { a.x += b.x; a.y += b.y; return a; }
inline Vec2f& operator-=(Vec2f& a, Vec2f b) { a.x -= b.x; a.y -= b.y; return a; }
inline Vec2f& operator*=(Vec2f& a, Vec2f b) { a.x *= b.x; a.y *= b.y; return a; }
inline Vec2f& operator/=(Vec2f& a, Vec2f b) { a.x /= b.x; a.y /= b.y; return a; }
inline Vec2f& operator*=(Vec2f& a, float s) { a.x *= s  ; a.y *= s  ; return a; }
inline Vec2f& operator/=(Vec2f& a, float s) { a.x /= s  ; a.y /= s  ; return a; }

inline Vec3f  operator+ (Vec3f  v)          { return { +v.x, +v.y, +v.z }; }
inline Vec3f  operator- (Vec3f  v)          { return { -v.x, -v.y, -v.z }; }
inline Vec3f  operator+ (Vec3f  a, Vec3f b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
inline Vec3f  operator- (Vec3f  a, Vec3f b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
inline Vec3f  operator* (Vec3f  a, Vec3f b) { return { a.x * b.x, a.y * b.y, a.z * b.z }; }
inline Vec3f  operator/ (Vec3f  a, Vec3f b) { return { a.x / b.x, a.y / b.y, a.z / b.z }; }
inline Vec3f  operator* (float  s, Vec3f a) { return { s   * a.x, s   * a.y, s   * a.z }; }
inline Vec3f  operator* (Vec3f  a, float s) { return { s   * a.x, s   * a.y, s   * a.z }; }
inline Vec3f  operator/ (float  s, Vec3f a) { return { s   / a.x, s   / a.y, s   / a.z }; }
inline Vec3f  operator/ (Vec3f  a, float s) { return { a.x / s  , a.y / s  , a.z / s   }; }
inline Vec3f& operator+=(Vec3f& a, Vec3f b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
inline Vec3f& operator-=(Vec3f& a, Vec3f b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; return a; }
inline Vec3f& operator*=(Vec3f& a, Vec3f b) { a.x *= b.x; a.y *= b.y; a.z *= b.z; return a; }
inline Vec3f& operator/=(Vec3f& a, Vec3f b) { a.x /= b.x; a.y /= b.y; a.z /= b.z; return a; }
inline Vec3f& operator*=(Vec3f& a, float s) { a.x *= s  ; a.y *= s  ; a.z *= s  ; return a; }
inline Vec3f& operator/=(Vec3f& a, float s) { a.x /= s  ; a.y /= s  ; a.z /= s  ; return a; }

inline float dot      (Vec2f a, Vec2f b) { return a.x * b.x + a.y * b.y; }
inline float len      (Vec2f v)          { return sqrtf(v.x * v.x + v.y * v.y); }
inline float sqrlen   (Vec2f v)          { return       v.x * v.x + v.y * v.y; }
inline Vec2f norm     (Vec2f v)          { return (1.0f / len(v)) * v; }
inline Vec2f norm_safe(Vec2f v)          { return (1.0f / (len(v) + 1e-8f)) * v; }
inline Vec2f abs      (Vec2f v)          { return { fabsf(v.x), fabsf(v.y) }; }
inline Vec2f cross    (Vec2f v)          { return { -v.y, v.x }; }

inline float dot      (Vec3f a, Vec3f b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline float len      (Vec3f v)          { return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z); }
inline float sqrlen   (Vec3f v)          { return       v.x * v.x + v.y * v.y + v.z * v.z; }
inline Vec3f norm     (Vec3f v)          { return (1.0f / len(v)) * v; }
inline Vec3f norm_safe(Vec3f v)          { return (1.0f / (len(v) + 1e-8f)) * v; }
inline Vec3f abs      (Vec3f v)          { return { fabsf(v.x), fabsf(v.y), fabsf(v.z) }; }
inline Vec3f cross    (Vec3f a, Vec3f b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }

template <typename T> inline void swap(T& a, T& b) { T t(a); a = b; b = t; }

namespace DOL
{
    inline Vec2f MinElts(const Vec2f& a, const Vec2f& b)
    {
        Vec2f result =
        {
            vl_min(a.x, b.x),
            vl_min(a.y, b.y)
        };
        return result;
    }

    inline Vec2f MaxElts(const Vec2f& a, const Vec2f& b)
    {
        Vec2f result =
        {
            vl_max(a.x, b.x),
            vl_max(a.y, b.y)
        };
        return result;
    }

    inline Vec3f MinElts(const Vec3f& a, const Vec3f& b)
    {
        Vec3f result =
        {
            vl_min(a.x, b.x),
            vl_min(a.y, b.y),
            vl_min(a.z, b.z)
        };
        return result;
    }

    inline Vec3f MaxElts(const Vec3f& a, const Vec3f& b)
    {
        Vec3f result =
        {
            vl_max(a.x, b.x),
            vl_max(a.y, b.y),
            vl_max(a.z, b.z)
        };
        return result;
    }

    inline int FloorToInt32(float x)
    {
        return int(floorf(x));
    }
    inline int RoundToInt32(float x)
    {
        return int(roundf(x));
    }
    inline int CeilToInt32(float x)
    {
        return int(ceilf(x));
    }

    struct RGBA8
    {
        union
        {
            struct { uint8_t r; uint8_t g; uint8_t b; uint8_t a; };
            uint8_t  c[4];
        };
    };

    struct Bounds2
    {
        Vec2f mMin;
        Vec2f mMax;
    };

    struct Bounds3
    {
        Vec3f mMin;
        Vec3f mMax;
    };

    inline Vec3f TriangleDoubleAreaNormal
    (
        const Vec3f& a,
        const Vec3f& b,
        const Vec3f& c
    )
    {
        Vec3f n =
        {
            (a.y - b.y) * (a.z + b.z) + (b.y - c.y) * (b.z + c.z) + (c.y - a.y) * (c.z + a.z),
            (a.z - b.z) * (a.x + b.x) + (b.z - c.z) * (b.x + c.x) + (c.z - a.z) * (c.x + a.x),
            (a.x - b.x) * (a.y + b.y) + (b.x - c.x) * (b.y + c.y) + (c.x - a.x) * (c.y + a.y)
        };

        return n;
    }

    inline Bounds2 TriangleBounds(const Vec2f& a, const Vec2f& b, const Vec2f& c)
    {
        Bounds2 box;

        box.mMin = MinElts(a, b);
        box.mMin = MinElts(box.mMin, c);

        box.mMax = MaxElts(a, b);
        box.mMax = MaxElts(box.mMax, c);

        return box;
    }

    inline Bounds3 TriangleBounds(const Vec3f& a, const Vec3f& b, const Vec3f& c)
    {
        Bounds3 box;

        box.mMin = MinElts(a, b);
        box.mMin = MinElts(box.mMin, c);

        box.mMax = MaxElts(a, b);
        box.mMax = MaxElts(box.mMax, c);

        return box;
    }
}

#endif
