Vec radiance(const Ray &ray, int depth, unsigned short *nRandSeeds) 
{
    double dist;                                   // distance to intersection
    int id = 0;                                 // id of intersected object
    if (!intersect(ray, dist, id)) return Vec();     // if miss, return black
    const Sphere &obj = spheres[id];            // the hit object
    Vec x = ray.o + ray.d * dist,                      // intersection point
        normalOut = (x - obj.p).norm(),                 // out-normal of sphere
        normal = normalOut.dot(ray.d) < 0 ? normalOut : normalOut * -1,       // fixed surface normal
        objColor = obj.c;
    double p = objColor.x > objColor.y && objColor.x > objColor.z ? objColor.x : objColor.y > objColor.z ? objColor.y : objColor.z; // max refl
    if (depth > 100)
        return obj.e; // Added to prevent stack overflow
    if (++depth > 5)
        if (erand48(nRandSeeds) < p)
            objColor = objColor * (1 / p);
        else
            return obj.e; //R.R.
    if (obj.refl == DIFF)           // Ideal DIFFUSE reflection
    {
        double r1 = 2 * pi * erand48(nRandSeeds), r2 = erand48(nRandSeeds), r2s = sqrt(r2);
        Vec w = normal, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + objColor.mult(radiance(Ray(x, d), depth, nRandSeeds));
    }
    else if (obj.refl == SPEC)      // Ideal SPECULAR reflection
    {
        return obj.e + objColor.mult(radiance(Ray(x, ray.d - normalOut * 2 * normalOut.dot(ray.d)), depth, nRandSeeds));
    }
    else
    {
        Ray reflRay(x, ray.d - normalOut * 2 * normalOut.dot(ray.d));           // Ideal dielectric REFRACTION
        bool into = normalOut.dot(normal) > 0;                          // Ray from outside going in?
        double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = ray.d.dot(normal), cos2t;
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)  // Total internal reflection
        {
            return obj.e + objColor.mult(radiance(reflRay, depth, nRandSeeds));
        }
        else
        {
            Vec tdir = (ray.d * nnt - normalOut * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
            double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(normalOut));
            double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
            // Russian roulette
            return obj.e + objColor.mult(depth > 2 ? (erand48(nRandSeeds) < P ? radiance(reflRay, depth, nRandSeeds) * RP
                                                               : radiance(Ray(x, tdir), depth, nRandSeeds) * TP)
                                            : radiance(reflRay, depth, nRandSeeds) * Re + radiance(Ray(x, tdir), depth, nRandSeeds) * Tr);
        }
    }
}