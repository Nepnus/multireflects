#include <stdint.h>
#include <math.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

EXPORT void core_thetas(const float* points_surface, const uint32_t* faces_ind, uint32_t faces_ind_len, float a, float b, float b2, float c,
    float* centroids, float* areas, float* unitNormals, float* thetas){
    
    const float PI = 3.141592653589793;
    const float PI_d2 = PI / 2.0;
    const float tol = 1.0e-6;

    uint32_t i, p1_ind, p2_ind, p3_ind;
    float tri_p1_x, tri_p1_y, tri_p1_z;
    float tri_p2_x, tri_p2_y, tri_p2_z;
    float tri_p3_x, tri_p3_y, tri_p3_z;
    float cen_x, cen_y, cen_z;
    float v1_x, v1_y, v1_z;
    float v2_x, v2_y, v2_z;
    float dirVect_x, dirVect_y, dirVect_z, dirVect_len;
    
    for(i=0; i<faces_ind_len; i++){
        p1_ind = faces_ind[i*3]; p2_ind = faces_ind[i*3+1]; p3_ind = faces_ind[i*3+2];

        tri_p1_x = points_surface[p1_ind*3]; tri_p1_y = points_surface[p1_ind*3+1] * b; tri_p1_z = points_surface[p1_ind*3+2] * b2; 
        tri_p2_x = points_surface[p2_ind*3]; tri_p2_y = points_surface[p2_ind*3+1] * b; tri_p2_z = points_surface[p2_ind*3+2] * b2; 
        tri_p3_x = points_surface[p3_ind*3]; tri_p3_y = points_surface[p3_ind*3+1] * b; tri_p3_z = points_surface[p3_ind*3+2] * b2;
        tri_p1_x = tri_p1_x < 0.0 ? tri_p1_x*a : tri_p1_x*c;
        tri_p2_x = tri_p2_x < 0.0 ? tri_p2_x*a : tri_p2_x*c;
        tri_p3_x = tri_p3_x < 0.0 ? tri_p3_x*a : tri_p3_x*c;

        cen_x = (tri_p1_x + tri_p2_x + tri_p3_x) / 3.0;
        cen_y = (tri_p1_y + tri_p2_y + tri_p3_y) / 3.0;
        cen_z = (tri_p1_z + tri_p2_z + tri_p3_z) / 3.0;
        centroids[i*3] = cen_x;
        centroids[i*3+1] = cen_y;
        centroids[i*3+2] = cen_z;

        if(cen_x < tol && cen_x > -tol)
            thetas[i] = PI_d2;
        else{
            v1_x = cen_y * cen_y;
            v1_y = cen_z * cen_z;
            v1_x += v1_y;
            v1_x = sqrtf(v1_x);
            v1_x /= cen_x;
            v1_x = atanf(v1_x);
            if(v1_x < 0.0)
                v1_x += PI;
            thetas[i] = v1_x;
        }

        v1_x = tri_p2_x - tri_p1_x;
        v1_y = tri_p2_y - tri_p1_y;
        v1_z = tri_p2_z - tri_p1_z;
        v2_x = tri_p3_x - tri_p1_x;
        v2_y = tri_p3_y - tri_p1_y;
        v2_z = tri_p3_z - tri_p1_z;

        dirVect_x = v1_y * v2_z - v1_z * v2_y;
        dirVect_y = v1_z * v2_x - v1_x * v2_z;
        dirVect_z = v1_x * v2_y - v1_y * v2_x;
        dirVect_len = sqrtf(dirVect_x*dirVect_x + dirVect_y*dirVect_y + dirVect_z*dirVect_z);
        dirVect_x /= dirVect_len; dirVect_y /= dirVect_len; dirVect_z /= dirVect_len;
        areas[i] = 0.5 * dirVect_len;

        dirVect_len = dirVect_x * cen_x + dirVect_y * cen_y + dirVect_z * cen_z;
        if(dirVect_len < 0.0){
            dirVect_x = - dirVect_x;
            dirVect_y = - dirVect_y;
            dirVect_z = - dirVect_z;
        }
        unitNormals[i*3] = dirVect_x;
        unitNormals[i*3+1] = dirVect_y;
        unitNormals[i*3+2] = dirVect_z;
    }
}

EXPORT void core_nothetas(const float* points_surface, const uint32_t* faces_ind, uint32_t faces_ind_len, float a, float b, float b2, float c,
    float* centroids, float* areas, float* unitNormals){
    
    uint32_t i, p1_ind, p2_ind, p3_ind;
    float tri_p1_x, tri_p1_y, tri_p1_z;
    float tri_p2_x, tri_p2_y, tri_p2_z;
    float tri_p3_x, tri_p3_y, tri_p3_z;
    float cen_x, cen_y, cen_z;
    float v1_x, v1_y, v1_z;
    float v2_x, v2_y, v2_z;
    float dirVect_x, dirVect_y, dirVect_z, dirVect_len;
    
    for(i=0; i<faces_ind_len; i++){
        p1_ind = faces_ind[i*3]; p2_ind = faces_ind[i*3+1]; p3_ind = faces_ind[i*3+2];

        tri_p1_x = points_surface[p1_ind*3]; tri_p1_y = points_surface[p1_ind*3+1] * b; tri_p1_z = points_surface[p1_ind*3+2] * b2; 
        tri_p2_x = points_surface[p2_ind*3]; tri_p2_y = points_surface[p2_ind*3+1] * b; tri_p2_z = points_surface[p2_ind*3+2] * b2; 
        tri_p3_x = points_surface[p3_ind*3]; tri_p3_y = points_surface[p3_ind*3+1] * b; tri_p3_z = points_surface[p3_ind*3+2] * b2;
        tri_p1_x = tri_p1_x < 0.0 ? tri_p1_x*a : tri_p1_x*c;
        tri_p2_x = tri_p2_x < 0.0 ? tri_p2_x*a : tri_p2_x*c;
        tri_p3_x = tri_p3_x < 0.0 ? tri_p3_x*a : tri_p3_x*c;

        cen_x = (tri_p1_x + tri_p2_x + tri_p3_x) / 3.0;
        cen_y = (tri_p1_y + tri_p2_y + tri_p3_y) / 3.0;
        cen_z = (tri_p1_z + tri_p2_z + tri_p3_z) / 3.0;
        centroids[i*3] = cen_x;
        centroids[i*3+1] = cen_y;
        centroids[i*3+2] = cen_z;

        v1_x = tri_p2_x - tri_p1_x;
        v1_y = tri_p2_y - tri_p1_y;
        v1_z = tri_p2_z - tri_p1_z;
        v2_x = tri_p3_x - tri_p1_x;
        v2_y = tri_p3_y - tri_p1_y;
        v2_z = tri_p3_z - tri_p1_z;

        dirVect_x = v1_y * v2_z - v1_z * v2_y;
        dirVect_y = v1_z * v2_x - v1_x * v2_z;
        dirVect_z = v1_x * v2_y - v1_y * v2_x;
        dirVect_len = sqrtf(dirVect_x*dirVect_x + dirVect_y*dirVect_y + dirVect_z*dirVect_z);
        dirVect_x /= dirVect_len; dirVect_y /= dirVect_len; dirVect_z /= dirVect_len;
        areas[i] = 0.5 * dirVect_len;

        dirVect_len = dirVect_x * cen_x + dirVect_y * cen_y + dirVect_z * cen_z;
        if(dirVect_len < 0.0){
            dirVect_x = - dirVect_x;
            dirVect_y = - dirVect_y;
            dirVect_z = - dirVect_z;
        }
        unitNormals[i*3] = dirVect_x;
        unitNormals[i*3+1] = dirVect_y;
        unitNormals[i*3+2] = dirVect_z;
    }
}


