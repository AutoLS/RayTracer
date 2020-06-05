#ifndef RAY_H
#define RAY_H

#define ArraySize(arr) (sizeof(arr)/sizeof(arr[0]))

enum light_type
{
	AMBIENT_LIGHT,
	DIRECTION_LIGHT,
	POINT_LIGHT,
};

struct light
{
	real32 Intensity;
	v3 Dir;
	v3 P;
	
	light_type Type;
};

struct material
{
	v3 Color;
	real32 Specular;
	real32 Reflective;
	bool Light;
};

struct sphere
{
	v3 P;
	real32 R;
	
	int MatIndex;
};

struct disk
{
	v3 P;
	v3 N;
	real32 r;
	
	int MatIndex;
};

struct intersection_result
{
	int HitIndex;
	real32 t;
	v3 P;
	v3 N;
};

struct world
{
	material* Materials;
	int MatCount;
	
	disk* Disks;
	int DiskCount;
	
	sphere* Spheres;
	int SphereCount;
	
	v3 ViewDir;
	
	light* Lights;
	int LightsCount;
};

#endif