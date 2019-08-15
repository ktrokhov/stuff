#include <iostream>
#include "stdio.h"
#include <cmath>
#include <cstdlib>

const int Cw = 800;
const int Ch = 800;
const int Vh = 1;
const int Vw = 1;
const double distance = 1;
const double T_MAX = 1000;
const double T_MIN = 1;
const char* FILE_NAME = "3.ppm";
const int RecursionDepth = 3;  
const int TransparencyDepth = 4; 
const int rotationX = 0;  
const int rotationY = 0;  
const int rotationZ = 0; 
const double Ox = 0;
const double Oy = 0;   
const double Oz = 0;   
const double Nref = 1; 
const int antialiasing_depth = 1; 




const double FogWidth = 8;  
int number = 0;

int scene = 0;
double FogRange = 5;
int shadow_depth = 20;




class Vector_ray; 
class Sphere;     
class Canvas;     
class Spheres;      
class Sphere_elem; 
class Light;       
class Lights;      



typedef struct Color {  
    int R, G, B;
    const Color operator*(double val) const {
        Color c;
        c.R = floor(double(R) * val);
        c.G = floor(double(G) * val);
        c.B = floor(double(B) * val);
        return c;
    }
    const Color operator+(Color C) const {
        Color K;
        K.R = R + C.R;
        K.G = G + C.G;
        K.B = B + C.B;
        return K;
    }
    const Color operator/(double val) const {
        Color c;
        c.R = floor(double(R) / val);
        c.G = floor(double(G) / val);
        c.B = floor(double(B) / val);
        return c;
    }
}Color;

const Color BACKGROUND_COLOR = {0, 184, 217}; 
const Color Fog = {133, 133, 133}; 

class Canvas {                 
    static Color mas[Ch * antialiasing_depth][Cw * antialiasing_depth];
public:
    Canvas() {
        for (int i = 0; i < Ch * antialiasing_depth; i++) {
            for (int j = 0; j < Cw * antialiasing_depth; j++) {
                mas[i][j].B = 255;
                mas[i][j].G = 255;
                mas[i][j].R = 255;
            }
        }
    }
    static void PutPixel(int x, int y, Color c) {
        mas[y - 1][x - 1].B = c.B;
        mas[y - 1][x - 1].G = c.G;
        mas[y - 1][x - 1].R = c.R;
    }
    static void Out() {
        FILE *f;
        f = fopen(FILE_NAME, "w+");
        fprintf(f, "P3\n%d %d\n255\n", Cw, Ch);
        for (int i = 0; i < Ch; i++) {
            for (int j = 0; j < Cw; j++) {
                fprintf(f, "%d %d %d\n", mas[i][j].R, mas[i][j].G, mas[i][j].B);
            }
        }
        fclose(f);
    }
    static void Antialiasing() { 
        Color average;
        average.R = 0;
        average.G = 0;
        average.B = 0;
        for (int p = 0; p < Ch * antialiasing_depth; p += antialiasing_depth) {
            for (int q = 0; q < Cw * antialiasing_depth; q += antialiasing_depth) {
                for (int k = p; k < p + antialiasing_depth; k++) {
                    for (int m = q; m < q + antialiasing_depth; m++) {
                        average = average + mas[k][m];
                    }
                }
                average.R = average.R / (antialiasing_depth * antialiasing_depth);
                average.G = average.G / (antialiasing_depth * antialiasing_depth);
                average.B = average.B / (antialiasing_depth * antialiasing_depth);
                mas[p / antialiasing_depth][q / antialiasing_depth] = average;
                average.R = 0;
                average.G = 0;
                average.B = 0;
            }    
        }
    }
};

class Vector_ray {         
    double x, y, z;
public:
    Vector_ray(double a, double b, double c): x(a), y(b), z(c) {}
    Vector_ray(): x(0), y(0), z(0) {};
    const Vector_ray operator+(const Vector_ray& b) {
        const Vector_ray k(x + b.x, y + b.y, z + b.z);
        return k;
    }
    const Vector_ray operator-(const Vector_ray& b) {
        const Vector_ray k(x - b.x, y - b.y, z - b.z);
        return k;
    }
    void SetPoints(double a, double b, double c) {
        x = a;
        y = b;
        z = c;
    }
    friend double dot(Vector_ray a, Vector_ray b) {
        return (a.x * b.x + a.y * b.y + a.z * b.z);
    }
    friend double length(Vector_ray a) {
        return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    }
    friend Vector_ray RotateCamera(Vector_ray D, double x, double y, double z);
    Vector_ray operator*(double val) {
        Vector_ray v;
        v.x = x * val;
        v.y = y * val;
        v.z = z * val;
        return v;
    }
    Vector_ray operator/(double val) {
        Vector_ray v;
        v.x = x / val;
        v.y = y / val;
        v.z = z / val;
        return v;
    }
    Vector_ray operator-() {
        Vector_ray v;
        v.x = -x;
        v.y = -y;
        v.z = -z;
        return v;
    }
    double ReturnX() {
        return x;
    }
    double ReturnY() {
        return y;
    }
    double ReturnZ() {
        return z;
    }
};

Color Canvas::mas[Ch * antialiasing_depth][Cw * antialiasing_depth];

class Sphere {         
    double radius;
    Vector_ray center;
    Color color_sphere;
    int specular;
    double reflective;
    double transparency;
    double refraction;
    double roughness;
    static int count;
public:
    Sphere(double x, double y, double z, double r, int Red, int Green, int Blue, 
        int specul, double reflect, double transp, double refract, double rough) {
        center.SetPoints(x, y, z);
        radius = r;
        color_sphere.R = Red;
        color_sphere.G = Green;
        color_sphere.B = Blue;
        specular = specul;
        reflective = reflect;
        transparency =transp;
        refraction = refract;
        roughness = rough;
        count++;
    }
    Sphere() {}
    double GetRadius() {
        return radius;
    }
    Vector_ray GetCenter() {
        return center;
    }
    Color GetSphereColor() {
        return color_sphere;
    }
    static int GetNumber() {
        return count;
    }
    int GetSpecular() {
        return specular;
    }
    double GetReflective() {
        return reflective;
    }
    double GetTransparensy() {
        return transparency;
    }
    double GetRefraction() {
        return refraction;
    }
    double GetRoughness() {
        return roughness;
    }
    ~Sphere() {
        count--;
    }
};

class Light {
    std::string type;
    double intensity;
    double radius;
    Vector_ray position_direction;
    static int count;
public:
    Light(std::string str, double inten, double x, double y, double z) {
        type = str;
        intensity = inten;
        position_direction.SetPoints(x, y, z);
        count++;
    }
    Light(std::string str, double inten) {
        type = str;
        intensity = inten;
        count++;
    }
    Light(std::string str, double inten, double x, double y, double z, double rad) {
        type = str;
        intensity = inten;
        position_direction.SetPoints(x, y, z);
        radius = rad;
        count++;
    }
    ~Light() {
        count--;
    }
    static int GetNumberLight() {
        return count;
    }
    std::string GetLightType() {
        return type;
    }
    double GetLightIntensity() {
        return intensity;
    }
    Vector_ray GetLightPosition() {
        return position_direction;
    }
    double GetLightRadius() {
        return radius;
    }
};

int Light::count = 0;

typedef struct Light_elem {
    Light_elem(std::string str, double inten, double x, double y, double z): elem(str, inten, x, y, z) {};
    Light_elem(std::string str, double inten): elem(str, inten) {};
    Light_elem(std::string str, double inten, double x, double y, double z, double radius): 
            elem(str, inten, x, y, z, radius) {};
    Light elem;
    Light_elem *next;
    
}Light_elem;

class Lights {
    static Light_elem* first;
    Light_elem *current;
public:
    void NewLight(std::string str, double inten, double x, double y, double z) {
        if (first == NULL) {
            first = new Light_elem(str, inten, x, y, z);
            first->next = NULL;
            current = first;
        } else {
            current->next = new Light_elem(str, inten, x, y, z);
            current = current->next;
            current->next = NULL;
        }
    }
    void NewLight(std::string str, double inten) {
        if (first == NULL) {
            first = new Light_elem(str, inten);
            first->next = NULL;
            current = first;
        } else {
            current->next = new Light_elem(str, inten);
            current = current->next;
            current->next = NULL;
        }
    }
    void NewLight(std::string str, double inten, double x, double y, double z, double rad) {
        if (first == NULL) {
            first = new Light_elem(str, inten, x, y, z, rad);
            first->next = NULL;
            current = first;
        } else {
            current->next = new Light_elem(str, inten, x, y, z, rad);
            current = current->next;
            current->next = NULL;
        }
    }
    ~Lights() {
        Light_elem *curr = first, *next_elem;
        if (first->next != NULL) {
            next_elem = first->next;
        } else {
            next_elem = NULL;
        }
        first = NULL;
        int k = Light::GetNumberLight();
        for (int i = 0; i < k; i++) {
            delete curr;
            if (next_elem != NULL) {
                curr = next_elem;
                next_elem = curr->next;
            }
        }
    }
    static std::string GetLightType(int k) {
        Light_elem *curr = first;
        for (int i = 0; i < k; i++) {
            curr = curr->next;
        }
        std::string str;
        str = curr->elem.GetLightType();
        return str;
    }
    static double GetLightIntensity(int k) {
        Light_elem *curr = first;
        for (int i = 0; i < k; i++) {
            curr = curr->next;
        }
        return curr->elem.GetLightIntensity();
    }
    static Vector_ray GetLightPosition(int k) {
        Light_elem *curr = first;
        for (int i = 0; i < k; i++) {
            curr = curr->next;
        }
        Vector_ray v = curr->elem.GetLightPosition();
        return v;
    }
    static Vector_ray GetLightDirection(int k) {
        Light_elem *curr = first;
        for (int i = 0; i < k; i++) {
            curr = curr->next;
        }
        Vector_ray v = curr->elem.GetLightPosition();
        return v;
    }
    static double GetLightRadius(int k) {
        Light_elem *curr = first;
        for (int i = 0; i < k; i++) {
            curr = curr->next;
        }
        double v = curr->elem.GetLightRadius();
        return v;
    }
};

Light_elem* Lights::first = NULL;

class Sphere_elem {          
public:
    Sphere elem;
    Sphere_elem* next;
    Sphere_elem(double x, double y, double z, double r, int Red, int Green, int Blue, int specul, 
        double reflect, double transp, double refract, double rough): 
    elem(Sphere(x, y, z, r, Red, Green, Blue, specul, reflect, transp, refract, rough)) {};  
};

int Sphere::count = 0;

class Spheres {                            
    static Sphere_elem *current;
    static Sphere_elem *first;
public:
    Spheres() {}
    static void NewSphere(double x, double y, double z, double r, int Red, int Green, int Blue, 
        int specul, double reflect, double transp, double refract, double rough) {
        if (first == NULL) {
            first = new Sphere_elem(x, y, z, r, Red, Green, Blue, specul, reflect, transp, refract, rough);
            current = first;
            current->next = NULL;
        } else {
            current->next = new Sphere_elem(x, y, z, r, Red, Green, Blue, specul, reflect, transp, refract, rough);
            current = current->next;
            current->next = NULL;
        }
    }
    static Color GetColor(int k) {
        Sphere_elem *curr = first;
        for (int i = 0; i < k - 1; i++) {
            curr = curr->next;
        }
        Color cl;
        cl = curr->elem.GetSphereColor();
        return cl;
    }
    static Vector_ray GetSphereCenter(int i) {
        Sphere_elem *curr = first;
        for (int k = 0; k < i - 1; k++) {
            curr = curr->next;
        }
        return curr->elem.GetCenter();
    }
    static double GetSphereRadius(int i) {
        Sphere_elem *curr = first;
        for (int k = 0; k < i - 1; k++) {
            curr = curr->next;
        }
        return curr->elem.GetRadius();
    }
    static int GetSpecular(int i) {
        Sphere_elem *curr = first;
        for (int k = 0; k < i - 1; k++) {
            curr = curr->next;
        }
        return curr->elem.GetSpecular();
    }
    static double GetReflective(int i) {
        Sphere_elem *curr = first;
        for (int k = 0; k < i - 1; k++) {
            curr = curr->next;
        }
        return curr->elem.GetReflective();
    }
    static double GetTransparensy(int i) {
        Sphere_elem *curr = first;
        for (int k = 0; k < i - 1; k++) {
            curr = curr->next;
        }
        return curr->elem.GetTransparensy();
    }
    static double GetRefraction(int i) {
        Sphere_elem *curr = first;
        for (int k = 0; k < i - 1; k++) {
            curr = curr->next;
        }
        return curr->elem.GetRefraction();
    }
    static double GetRoughness(int i) {
        Sphere_elem *curr = first;
        for (int k = 0; k < i - 1; k++) {
            curr = curr->next;
        }
        return curr->elem.GetRoughness();
    }
    ~Spheres() {
        Sphere_elem *curr = first, *next_elem;
        if (first->next != NULL) {
            next_elem = first->next;
        } else {
            next_elem = NULL;
        }
        first = NULL;
        int k = Sphere::GetNumber();
        for (int i = 0; i < k; i++) {
            delete curr;
            if (next_elem != NULL) {
                curr = next_elem;
                next_elem = curr->next;
            }
        }
    }
};



Sphere_elem* Spheres::first = NULL;
Sphere_elem* Spheres::current = NULL;

void IntersectRaySphere(Vector_ray O, Vector_ray D, int i, double &t1, double &t2) { //_____Решение квадратного уравнения____
    Vector_ray C, oc;
    double k1, k2, k3, discriminant, r;
    C = Spheres::GetSphereCenter(i); 
    r = Spheres::GetSphereRadius(i);
    oc = O - C;
    k1 = dot(D, D);
    k2 = 2 * dot(oc, D);
    k3 = dot(oc, oc) - r * r;
    discriminant = k2 * k2 - 4 * k1 * k3;
    if (discriminant < 0) {
        t1 = T_MAX;
        t2 = T_MAX;
    } else {
        t1 = (-k2 + sqrt(discriminant)) / (2 * k1);
        t2 = (-k2 - sqrt(discriminant)) / (2 * k1);
    }
}

void ClosestIntersection(Vector_ray O, Vector_ray D, double t_min, double t_max, int& closest_sphere, double& closest_t) {
    double t1, t2;
    int k = Sphere::GetNumber();
    int closest_sph = closest_sphere;
    for (int i = 0; i < k; i++) {
        IntersectRaySphere(O, D, i + 1, t1, t2);
        if (((t1 < t_max) && (t1 > t_min)) && (t1 < closest_t)) {
            closest_t = t1;
            closest_sphere = i + 1;
        }
        if (((t2 < t_max) && (t2 > t_min)) && (t2 < closest_t)) {
            closest_t = t2;
            closest_sphere = i + 1;
        }
    }
}

Vector_ray ReflectRay(Vector_ray R, Vector_ray N) {  
    return (N * 2 * dot(N, R) - R);
}

Vector_ray RotateVector(Vector_ray V, Vector_ray L, double phi) {  
    Vector_ray result;
    L = L / length(L);
    double resX, resY, resZ;
    double Vx = V.ReturnX();
    double Vy = V.ReturnY();
    double Vz = V.ReturnZ();
    double Lx = L.ReturnX();
    double Ly = L.ReturnY();
    double Lz = L.ReturnZ();
    resX = Vx * (cos(phi) + (1 - cos(phi)) * Lx * Lx) + Vy * ((1 - cos(phi)) * Ly * Lx + sin(phi) * Lz) + Vz * ((1 - cos(phi)) * Lz * Lx - sin(phi) * Ly);
    resY = Vx * ((1 - cos(phi)) * Lx * Ly - sin(phi) * Lz) + Vy * (cos(phi) + (1 - cos(phi)) * Ly * Ly) + Vz * ((1 - cos(phi)) * Lz * Ly + sin(phi) * Lx);
    resZ = Vx * ((1 - cos(phi)) * Lz * Lx + sin(phi) * Ly) + Vy * ((1 - cos(phi)) * Ly * Lz - sin(phi) * Lx) + Vz * (cos(phi) + (1 - cos(phi)) * Lz * Lz);
    result.SetPoints(resX, resY, resZ);
    return result;
}

Vector_ray ComputeVectorRadius(Vector_ray L, Vector_ray center, double radius) { 
    double d; 
    Vector_ray A, OA, Rad;
    Vector_ray O_X, O_Y, O_Z;
    O_X.SetPoints(1, 0, 0);
    O_Y.SetPoints(0, 1, 0);
    O_Z.SetPoints(0, 0, 1);
    d = -dot(L, center);
    if ((dot(L, O_X) == 0) && (dot(L, O_Z) == 0)) {
        A.SetPoints(0, -d / L.ReturnY(), 0);
    } else if ((dot(L, O_Y) == 0) && (dot(L, O_Z) == 0)) {
        A.SetPoints(-d / L.ReturnX(), 0, 0);
    } else {
        A.SetPoints(0, 0, -d / L.ReturnZ());
    }
    OA = A - center;
    OA = OA / length(OA);
    Rad = OA * radius * double(rand()) / double(RAND_MAX);
    double phi = 2 * M_PI * double(rand()) / double(RAND_MAX);
    Rad = RotateVector(Rad, L, phi);
    return Rad;
}

double ComputePointLighting(Vector_ray L, Vector_ray center, Vector_ray P, Vector_ray N, Vector_ray V, int closest_sphere, double radius, double light_intens) { //____Освещенность от точечного источника____
    Vector_ray R, L2; 
    double intens = 0;
    double N_dot_L, R_dot_V;
    for (int i = 0; i < shadow_depth; i++) {
        L2 = L + ComputeVectorRadius(L, center, radius);
        int shadow_sphere = -1;   
        double shadow_t = 1;
        ClosestIntersection(P, L2, 0.001, 1, shadow_sphere, shadow_t);
        if (shadow_sphere != -1) {
            continue;
        }
        N_dot_L = dot(N, L2); 
        if (N_dot_L > 0) {
            intens += light_intens * N_dot_L / (length(N) * length(L2));
        }
        if (Spheres::GetSpecular(closest_sphere) != -1) { 
            R = ReflectRay(L2, N);
            R_dot_V = dot(R, V);
            if (R_dot_V > 0) {
                intens += light_intens * pow(R_dot_V / (length(R) * length(V)), Spheres::GetSpecular(closest_sphere));
            }
        }
    }
    intens = intens / shadow_depth;
    return intens;
}

double ComputeLighting(Vector_ray P, Vector_ray N, Vector_ray V, int closest_sphere) {  
    double intens = 0, N_dot_L, R_dot_V, t_max, point_intens;
    Vector_ray L, R;
    int k = Light::GetNumberLight();
    for (int j = 0; j < k; j++) {
        if (Lights::GetLightType(j) == "ambient") {
            intens += Lights::GetLightIntensity(j);
        } else {
            if (Lights::GetLightType(j) == "point") {
                L = Lights::GetLightPosition(j) - P;
                //t_max = 1;
                point_intens = ComputePointLighting(L, Lights::GetLightPosition(j), P, N, V, closest_sphere, Lights::GetLightRadius(j), Lights::GetLightIntensity(j));
                intens += point_intens;
                continue;
            } else {
                L = Lights::GetLightDirection(j);
                t_max = T_MAX;
            }
            int shadow_sphere = -1; 
            double shadow_t = t_max;
            ClosestIntersection(P, L, 0.001, t_max, shadow_sphere, shadow_t);
            if (shadow_sphere != -1) {
                continue;
            }
            N_dot_L = dot(N, L);  
            if (N_dot_L > 0) {
                intens += Lights::GetLightIntensity(j) * N_dot_L / (length(N) * length(L));
            }
            if (Spheres::GetSpecular(closest_sphere) != -1) { 
                R = ReflectRay(L, N);
                R_dot_V = dot(R, V);
                if (R_dot_V > 0) {
                    intens += Lights::GetLightIntensity(j) * pow(R_dot_V / (length(R) * length(V)), Spheres::GetSpecular(closest_sphere));
                }
            }
        }
    }
    return intens;
}

Vector_ray IntersectRayRefract(Vector_ray D, Vector_ray N, double n1, double n2) {
    double p1, p2, p3;
    double N_dot_D = dot(N, D);
    double N_dot_N = dot(N, N);
    double D_dot_D = dot(D, D);
    double t1, t2;
    if ((n1 == n2) || (N_dot_D / (length(D) * length(N)) == 1)) {
        return D;
    }
    p1 = (1 - n1 * n1 / (n2 * n2)) * (N_dot_D * N_dot_D - N_dot_N * D_dot_D);
    p2 = 2 * n1 * n1 / (n2 * n2) * (N_dot_D * N_dot_N - N_dot_D * N_dot_D * N_dot_D / D_dot_D);
    p3 = n1 * n1 / (n2 * n2) * (N_dot_N * N_dot_N - N_dot_N * N_dot_D * N_dot_D / D_dot_D);
    double discriminant = p2 * p2 - 4 * p1 * p3;
    if (discriminant < 0) {
        return D;
    }
    t1 = (-p2 + sqrt(discriminant)) / (2 * p1);
    t2 = (-p2 - sqrt(discriminant)) / (2 * p1);
    if (t1 > 0) {
        return N + D * t1;
    } else if (t2 > 0) {
        return N + D * t2;
    } else {
        return D;
    }
}

Color TraceRay(Vector_ray O, Vector_ray D, double t_min, double t_max, int recursion_depth, int transparency_depth, int closest_sph) {  
    int closest_sphere = -1;
    double closest_t = t_max;
    Color local_color;
    Vector_ray R;
    ClosestIntersection(O, D, t_min, t_max, closest_sphere, closest_t);
    if (closest_sphere == closest_sph) {
        Vector_ray P = O + D * closest_t;
        Vector_ray N = P - Spheres::GetSphereCenter(closest_sphere);
        N = N / length(N);
        D = IntersectRayRefract(D, N, Spheres::GetRefraction(closest_sphere), Nref);
        t_min = closest_t;
        closest_sphere = -1;
        closest_t = t_max;
        ClosestIntersection(O, D, t_min, t_max, closest_sphere, closest_t);
    }
    if (closest_sphere == -1) {
        return Fog;
        
    } else {
        Color result_color;
        Vector_ray P = O + D * closest_t;
        if (P.ReturnZ() > FogRange + FogWidth) {
            return Fog;
        }
        Vector_ray N = P - Spheres::GetSphereCenter(closest_sphere);
        N = N / length(N);                            
        double Radius = Spheres::GetRoughness(closest_sphere);
        Vector_ray Q = P + N;                      
        N = N + ComputeVectorRadius(N, Q, Radius); 
        if(number == 1) {
            if (closest_sphere == 1) {   
                if (((int(abs(P.ReturnX())) % 2 + int(abs(P.ReturnZ())) % 2) % 2 == 1) && (P.ReturnX() >= 0)) {
                    closest_sphere = 2;
                } else if (((int(abs(P.ReturnX())) % 2 + int(abs(P.ReturnZ())) % 2) % 2 == 0) && (P.ReturnX() < 0)) {
                    closest_sphere = 2;
                }
            }
        }
        if (ComputeLighting(P, N, -D, closest_sphere) > 1) {  
            local_color = Spheres::GetColor(closest_sphere);
        } else {
            local_color = (Spheres::GetColor(closest_sphere)) * ComputeLighting(P, N, -D, closest_sphere);
        }
        double refl;                                   
        refl = Spheres::GetReflective(closest_sphere);
        double transparency = Spheres::GetTransparensy(closest_sphere);
        
        R = ReflectRay(-D, N); 
        if ((recursion_depth > 0) && (refl > 0)) { 
            Color reflected_color = TraceRay(P, /*RandRefl*/R, 0.001, T_MAX, recursion_depth - 1, transparency_depth, closest_sphere);
            result_color = local_color * (1 - refl) + reflected_color * refl;
        } else {
            result_color = local_color;
        }
        if ((transparency_depth > 0) && (transparency > 0)) { 
            D = IntersectRayRefract(D, -N, Nref, Spheres::GetRefraction(closest_sphere));
            Color transparency_color = TraceRay(P, D, 0.00001, T_MAX, recursion_depth, transparency_depth - 1, closest_sphere);
            result_color = result_color * (1 - transparency) + transparency_color * transparency;
        }
        if (P.ReturnZ() > FogRange) {                    
            double k = (P.ReturnZ() - FogRange) / FogWidth;
            if (k > 1) {
                k = 1;
            }
            result_color = Fog * k + result_color * (1 - k);
        }
        return result_color;    
    }
}

Vector_ray CanvasToViewport(int x, int y) {   
    Vector_ray k(double(x) * Vw / (Cw * antialiasing_depth), double(y) * Vh / (Ch * antialiasing_depth), distance);
    return k;
}

Vector_ray RotateCamera(Vector_ray D, double angle_x, double angle_y, double angle_z) {
    D.SetPoints(D.x, D.y * cos(angle_x) + D.z * sin(angle_x), D.z * cos(angle_x) - D.y * sin(angle_x)); 
    D.SetPoints(D.x * cos(angle_y) - D.z * sin(angle_y), D.y, D.x * sin(angle_y) + D.z * cos(angle_y)); 
    D.SetPoints(D.x * cos(angle_z) + D.y * sin(angle_z), D.y * cos(angle_z) - D.x * sin(angle_z), D.z); 
    return D;   
}

void GetImage() {                    
    Vector_ray O, D;
    O.SetPoints(Ox, Oy, Oz);
    Color cl;
    int process1 = 0;
    int i = 0;
    for (int x = -(Cw * antialiasing_depth) / 2; x <= (Cw * antialiasing_depth) / 2; x++) {
        for (int y = -(Ch * antialiasing_depth) / 2; y <= (Ch * antialiasing_depth) / 2; y++) {
            D = CanvasToViewport(x, y);
            D = RotateCamera(D, rotationX * M_PI / 180, rotationY * M_PI / 180, rotationZ * M_PI / 180);
            if ((((Cw * antialiasing_depth) / 2 + x > 0) && 
               ((Cw * antialiasing_depth) / 2 + x <= (Cw * antialiasing_depth))) &&
               (((Ch * antialiasing_depth) / 2 - y > 0) && 
               ((Ch * antialiasing_depth) / 2 - y <= (Ch * antialiasing_depth)))) {
                cl = TraceRay(O, D, 1, T_MAX, RecursionDepth, TransparencyDepth, -2);
                Canvas::PutPixel((Cw * antialiasing_depth) / 2 + x, (Ch * antialiasing_depth) / 2 - y, cl);
            }
        }
        int process = ceil(100 * i / double(Ch * antialiasing_depth)) / 10;  
        if (process != process1) {
           std::cout << "Идет процесс " << process * 10 << "% " << std::endl;
        }
        process1 = process;
        i++;
    }
    Canvas::Antialiasing();
}






int main() {

    std::cout << "введите номер сцены:\n";
    std::cin >> scene;
    switch(scene) {
        case 1:
        {
            Canvas c;
            Spheres s;   
            Lights l;    
            FogRange = 11.9;
            shadow_depth = 15;
            
            Spheres::NewSphere(0, -1, 6, 1, 255, 0, 0, 500, 0.2, 0, Nref,0); 
            Spheres::NewSphere(2, 0, 7, 1, 0, 0, 255, 500, 0.3, 0, Nref,0);
            Spheres::NewSphere(-2, 0, 7, 1, 0, 255, 0, 30, 1, 0, Nref,0); 
            Spheres::NewSphere(2, 1, 9, 2, 255, 255, 0, 30, 0.2, 0, Nref,0);
            Spheres::NewSphere(-2, 2, 9, 1.3, 128, 0, 128, 500, 0.1, 0.2, Nref,0); 
            Spheres::NewSphere(0, -5001, 0, 5000, 255, 255, 0, 10, 0, 0, Nref,0); 
            Spheres::NewSphere(0, 0, 5012, 5000, 0, 200, 0, 1, 0, 0, Nref,0); 
            Spheres::NewSphere(-5004, 0, 0, 5000, 0, 200, 0, -1, 0, 0, Nref,0); 
            Spheres::NewSphere(5004, 0, 0, 5000, 0, 200, 0, -1, 0, 0, Nref,0); 
            Spheres::NewSphere(0, 5004, 0, 5000, 0, 200, 0, -1, 0, 0, Nref,0); 
            Spheres::NewSphere(0, 0, -5001, 5000, 0, 200, 0, -1, 0, 0, Nref,0); 
            
            
            l.NewLight("ambient", 0.1);
            l.NewLight("point", 0.8, 2, 1, 2, 0.5);
            l.NewLight("directional", 0.1, 1, 4, 4);
            
            GetImage();
            Canvas::Out();

            break;
        }
        case 2:
        {
            Canvas c;
            Spheres s;   
            Lights l;    

            
            
            FogRange = 10;
            number = 1;
            Spheres::NewSphere(0, -10001.5, 8, 10000, 255, 255, 255, -1, 0, 0, Nref, 0.05); 
            Spheres::NewSphere(0, -10001.5, 8, 10000, 0, 0, 0, -1, 0, 0, Nref, 0); 
            Spheres::NewSphere(-1, -0.5, 5, 1, 255, 0, 0, 30, 0.1, 0.6, 1.33, 0); 
            Spheres::NewSphere(0.5, -1, 4, 0.5, 0, 255, 0, 500, 0.2, 0.8, 1.33, 0); 
            Spheres::NewSphere(0.5, 0, 9.5, 1.5, 0, 0, 255, 30, 0, 0.7, 1.33, 0); 

            
            l.NewLight("ambient", 0.1);
            l.NewLight("point", 0.8, 2, 1, 2, 0.5);
            l.NewLight("directional", 0.1, 1, 4, 4);
            
            GetImage();
            Canvas::Out();

            
            break;
        }
        case 3:
        {
            Canvas c;
            Spheres s;  
            Lights l;   

            
            
            shadow_depth = 10;
            Spheres::NewSphere(0, -10001.5, 8, 10000, 255, 255, 0, -1, 0, 0, Nref, 0.05); 
            Spheres::NewSphere(0.5, -1, 4, 1, 255, 0, 0, 30, 0.4, 0, Nref, 0.01); 
            Spheres::NewSphere(-1, -0.5, 5, 0.5, 0, 255, 0, 500, 0.2, 0.8, 1.33, 0); 
            Spheres::NewSphere(-3, 0, 2, 1.5, 0, 0, 255, 10, 0, 0, Nref, 0); 

            l.NewLight("ambient", 0.1);
            l.NewLight("point", 0.8, 2, 1, 2, 0);
            l.NewLight("directional", 0.1, 1, 4, 4);
            
            GetImage();
            Canvas::Out();

            
            break;
        }
        default:
        {
            
            std::cout << "неверный ввод\n";
            break;
        }
    }
    
    return 0;
}
