#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

void Intersection::computePartials(const RayDifferential &ray) {
    Float A[2][2], Bx[2], By[2], x[2];
    int axes[2];

    /* Compute the texture coordinates partials wrt.
       changes in the screen-space position. Based on PBRT */
    if (hasUVPartials || !ray.hasDifferentials)
        return;

    hasUVPartials = true;

    if (dpdu.isZero() && dpdv.isZero()) {
        dudx = dvdx = dudy = dvdy = 0.0f;
        return;
    }

    /* Compute a few projections onto the surface normal */
    const Float
        pp  = dot(geoFrame.n, Vector(p)),
        pox = dot(geoFrame.n, Vector(ray.rxOrigin)),
        poy = dot(geoFrame.n, Vector(ray.ryOrigin)),
        prx = dot(geoFrame.n, ray.rxDirection),
        pry = dot(geoFrame.n, ray.ryDirection);

    if (EXPECT_NOT_TAKEN(prx == 0 || pry == 0)) {
        dudx = dvdx = dudy = dvdy = 0.0f;
        return;
    }

    /* Compute ray-plane intersections against the offset rays */
    const Float tx = (pp-pox) / prx, ty = (pp-poy) / pry;

    /* Calculate the U and V partials by solving two out
       of a set of 3 equations in an overconstrained system */
    Float absX = std::abs(geoFrame.n.x),
          absY = std::abs(geoFrame.n.y),
          absZ = std::abs(geoFrame.n.z);

    if (absX > absY && absX > absZ) {
        axes[0] = 1; axes[1] = 2;
    } else if (absY > absZ) {
        axes[0] = 0; axes[1] = 2;
    } else {
        axes[0] = 0; axes[1] = 1;
    }

    A[0][0] = dpdu[axes[0]];
    A[0][1] = dpdv[axes[0]];
    A[1][0] = dpdu[axes[1]];
    A[1][1] = dpdv[axes[1]];

    /* Auxilary intersection point of the adjacent rays */
    Point px = ray.rxOrigin + ray.rxDirection * tx,
          py = ray.ryOrigin + ray.ryDirection * ty;

    Bx[0] = px[axes[0]] - p[axes[0]];
    Bx[1] = px[axes[1]] - p[axes[1]];
    By[0] = py[axes[0]] - p[axes[0]];
    By[1] = py[axes[1]] - p[axes[1]];

    if (EXPECT_TAKEN(solveLinearSystem2x2(A, Bx, x))) {
        dudx = x[0]; dvdx = x[1];
    } else {
        dudx = 1; dvdx = 0;
    }

    if (EXPECT_TAKEN(solveLinearSystem2x2(A, By, x))) {
        dudy = x[0]; dvdy = x[1];
    } else {
        dudy = 0; dudy = 1;
    }
}

std::string Intersection::toString() const {
    if (!isValid())
        return "Intersection[invalid]";
    std::ostringstream oss;
    oss << "Intersection[" << endl
        << "  p = " << p.toString() << "," << endl
        << "  wi = " << wi.toString() << "," << endl
        << "  t = " << t << "," << endl
        << "  geoFrame = " << indent(geoFrame.toString()) << "," << endl
        << "  shFrame = " << indent(shFrame.toString()) << "," << endl
        << "  uv = " << uv.toString() << "," << endl
        << "  hasUVPartials = " << hasUVPartials << "," << endl
        << "  dpdu = " << dpdu.toString() << "," << endl
        << "  dpdv = " << dpdv.toString() << "," << endl;
    if (hasUVPartials) {
        oss << "  dud[x,y] = [" << dudx << ", " << dudy << "]," << endl
            << "  dvd[x,y] = [" << dvdx << ", " << dvdy << "]," << endl;
    }
    oss << "  time = " << time << "," << endl
        << "  shape = " << indent(((Object *)shape)->toString()) << endl
        << "]";
    return oss.str();
}

MTS_NAMESPACE_END
