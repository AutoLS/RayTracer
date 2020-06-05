#include <stdlib.h>
#include <stdio.h>
#include <direct.h>

#include "AE/math.h"
#include "ray.h"

#pragma pack(push, 1)
struct bitmap_header
{
	uint16 FileType;
	uint32 FileSize;
	uint16 Reserved1;
	uint16 Reserved2;
	uint32 BitmapOffset;
	uint32 Size;
	int Width;
	int Height;
	uint16 Planes;
	uint16 BitsPerPixel;
	uint32 Compression;
	uint32 SizeOfBitmap;
	int HorResolution;
	int VertResolution;
	uint32 ColorsUsed;
	uint32 ColorsImportant;	
};
#pragma pack(pop)

struct image
{
	int Width;
	int Height;
	uint32* Pixels;
};

void WriteImage(char* FileName, image* Image)
{
	int OutputPixelSize = Image->Width*Image->Height*sizeof(int);
	
	bitmap_header Header = {};
	Header.FileType = 0x4D42;
	Header.FileSize = sizeof(Header) + OutputPixelSize;
	Header.BitmapOffset = sizeof(Header);
	Header.Size = sizeof(Header) - 14;
	Header.Width = Image->Width;
	Header.Height = Image->Height;
	Header.Planes = 1;
	Header.BitsPerPixel = 32;
	Header.Compression = 0;
	Header.SizeOfBitmap = OutputPixelSize;
	Header.HorResolution = 0;
	Header.VertResolution = 0;
	
	char Path[100] = {};
	sprintf(Path, "output/%s.bmp", FileName);
	
	FILE* OutFile = fopen(Path, "wb");
	if(OutFile)
	{
		fwrite(&Header, sizeof(Header), 1, OutFile);
		fwrite(Image->Pixels, OutputPixelSize, 1, OutFile);
		fclose(OutFile);
		printf("Writting image successfully...\n");
	}
	else
	{
		fprintf(stderr, "[ERROR] Unable to write output file %s.\n", Path);
	}
}

void FillImage(image* Image, int R, int G, int B)
{
	uint32* Pixels = Image->Pixels;
	for(int y = 0; y < Image->Height; ++y)
	{
		for(int x = 0; x < Image->Width; ++x)
		{
			*Pixels++ = R << 16 | G << 8 | B << 0;
		}
	}
}

intersection_result 
RayIntersection(world* World, 
					v3 RayOrigin, v3 RayDir, 
					real32 t_min, real32 t_max)
{
	intersection_result Result = {};
	Result.t = FLT_MAX;
	Result.HitIndex = -1;
	
	for(int i = 0; i < World->SphereCount; ++i)
	{
		sphere Sphere = World->Spheres[i];
		
		v3 OC = RayOrigin - Sphere.P;
		
		real32 a = Dot(RayDir, RayDir);
		real32 b = 2 * Dot(OC, RayDir);
		real32 c = Dot(OC, OC) - (Sphere.R*Sphere.R);
		
		real32 RootTerm = sqrtf((b*b) - 4*a*c);
		real32 Denom = 2 * a;
		if(RootTerm > 0)
		{
			real32 t1 = (-b + RootTerm) / Denom;
			real32 t2 = (-b - RootTerm) / Denom;
			if(t2 < t1) t1 = t2;
			if((t1 >= t_min && t1 <= t_max) && t1 < Result.t)
			{
				Result.t = t1;
				Result.HitIndex = Sphere.MatIndex;
				
				Result.P = RayOrigin + Result.t*RayDir;
				Result.N = Normalize(Result.P - Sphere.P);
			}
		}
	}
	
	for(int i = 0; i < World->DiskCount; ++i)
	{
		disk Disk = World->Disks[i];
		
		real32 Denom = Dot(Disk.N, RayDir);
		if(Denom > 0)
		{
			real32 t = Dot(Disk.P - RayOrigin, Disk.N) / Denom;
			if((t >= t_min && t <= t_max) && t < Result.t)
			{
				v3 P = RayOrigin + t*RayDir;
				real32 Distance = Length(P - Disk.P);
				if(Distance <= Disk.r)
				{
					Result.t = t;
					Result.HitIndex = Disk.MatIndex;
					
					Result.P = P;
					Result.N = Disk.N;
				}
			}
		}
	}
	
	return Result;
}

real32 ComputeLight(world* World, v3 P, v3 N, real32 Shininess)
{
	real32 Result = 0;
	
	v3 L = {};
	real32 t_max = FLT_MAX;
	
	for(int i = 0; i < World->LightsCount; ++i)
	{
		light Light = World->Lights[i];
		if(Light.Type == AMBIENT_LIGHT)
		{
			Result += Light.Intensity;
		}
		else
		{
			switch(Light.Type)
			{
				case DIRECTION_LIGHT:
				{
					L = Normalize(Light.Dir);
				} break;
				case POINT_LIGHT:
				{
					L = Normalize(Light.P - P);
					t_max = 1.0f;
				} break;
			}
			
			intersection_result Shadow = 
			RayIntersection(World, P, L, 0.001f, t_max);
			if(Shadow.HitIndex >= 0 && !World->Materials[Shadow.HitIndex].Light) 
				continue;
			
			real32 Diffuse = Dot(L, N);
			if(Diffuse > 0)
			Result += Light.Intensity * Diffuse;
		
			v3 R = Normalize(Reflect(-L, N));
			real32 Specular = Dot(R, World->ViewDir);
			if(Specular > 0)
			{
				real32 Reflection = (real32)pow(Specular, Shininess);
				Result += Light.Intensity * Reflection;
			}
		}
	}

	return Result;
}

v3 RayTrace(v3 RayOrigin, v3 RayDir, world* World, 
			real32 t_min, real32 t_max, int Depth)
{
	v3 LocalColor = {}, ReflectColor = {};
	real32 r;
	World->ViewDir = Normalize(-RayDir);
	
	intersection_result Intersection = 
	RayIntersection(World, RayOrigin, RayDir, t_min, t_max);
	
	int Index = Intersection.HitIndex;
	
	if(Index >= 0)
	{
		material Mat = World->Materials[Index];
		v3 P = Intersection.P;
		v3 N = Intersection.N;
				
		real32 Intensity = 1.0f;
		if(!Mat.Light)
		{
			Intensity = ComputeLight(World, P, N, Mat.Specular);
			if(Intensity > 1.0f)
				Intensity = 1.0f;
		}
		
		LocalColor = Mat.Color * Intensity;
		
		r = Mat.Reflective;
		if(Depth <= 0 || r <= 0)
		{
			return LocalColor;
		}
		
		v3 R = Reflect(RayDir, N);
		ReflectColor = RayTrace(P, R, World, 0.001f, FLT_MAX, Depth - 1);
	}
	
	return Lerp(LocalColor, ReflectColor, r);
}

static real32
ExactLinearTosRGB(real32 L)
{
	if(L < 0)
		L = 0;
	if(L > 1.0f)
		L = 1.0f;
	
	real32 S = L * 12.92f;
	if(L > 0.0031308)
	{
		S = 1.055f * powf(L, 1.0f/2.4f) - 0.055f;
	}
	
	return S;
}

int main(int ArgCount, char** Args)
{	
	if(ArgCount == 2)
	{
		_mkdir("output");
		
		image Image = {};
		Image.Width = 1920;//atoi(Args[1]);
		Image.Height = 1080;//atoi(Args[2]);
		printf("Width: %d Height: %d\n", Image.Width, Image.Height);
		Image.Pixels = (uint32*)malloc(sizeof(int) * Image.Width * Image.Height);
		
		material Materials[6] = {};
		Materials[0].Color = V3(1, 0, 0);
		Materials[0].Specular = 500;
		Materials[0].Reflective = 0;
		
		Materials[1].Color = V3(0.6f, 1.0f, 0.3f);
		Materials[1].Specular = 10;
		Materials[1].Reflective = 0.2f;
		
		Materials[2].Color = V3(0.1f, 0.8f, 1);
		Materials[2].Specular = 500;
		Materials[2].Reflective = 0.7f;
		
		Materials[3].Color = V3(0.5f, 0.5f, 0.5f);
		Materials[3].Specular = 1000;
		Materials[3].Reflective = 0;
		
		Materials[4].Color = V3(1, 1, 1);
		Materials[4].Specular = 0;
		Materials[4].Reflective = 0;
		Materials[4].Light = true;
		
		Materials[5].Color = V3(1, 1, 1);
		Materials[5].Specular = 1000;
		Materials[5].Reflective = 1.0f;
		
		sphere Spheres[5] = {};
		Spheres[0].P = V3(0, -1, 3);
		Spheres[0].R = 1.0f;
		Spheres[0].MatIndex = 0;
		
		Spheres[1].P = V3(-2, 0, 4);
		Spheres[1].R = 1;
		Spheres[1].MatIndex = 1;
		
		Spheres[2].P = V3(2, 0, 4);
		Spheres[2].R = 1;
		Spheres[2].MatIndex = 2;
		
		Spheres[3].P = V3(0, -5001, 0);
		Spheres[3].R = 5000;
		Spheres[3].MatIndex = 3;
		
		Spheres[4].P = V3(0, 1, 0);
		Spheres[4].R = 0.5f;
		Spheres[4].MatIndex = 4;
		
		disk Disks[1] = {};
		Disks[0].P = V3(0, 2, 4);
		Disks[0].N = Normalize(Disks[0].P);
		Disks[0].r = 1.5f;
		Disks[0].MatIndex = 5;
		
		light Lights[3] = {};
		Lights[0].Intensity = 0.2f;
		Lights[0].Type = AMBIENT_LIGHT;
		
		Lights[1].Intensity = 0.0f;
		Lights[1].Dir = V3(1, 4, 4);
		Lights[1].Type = DIRECTION_LIGHT;
		
		Lights[2].Intensity = 0.6f;
		Lights[2].P = Spheres[4].P;
		Lights[2].Type = POINT_LIGHT;
		
		world World = {};
		World.Materials = Materials;
		World.MatCount = ArraySize(Materials);
		World.Spheres = Spheres;
		World.SphereCount = ArraySize(Spheres);
		World.Disks = Disks;
		World.DiskCount = ArraySize(Disks);
		World.ViewDir = V3(0, 0, -1);
		World.Lights = Lights;
		World.LightsCount = ArraySize(Lights);
		
		v2 FilmDim = {2.0f, 2.0f};
		if(Image.Width > Image.Height)
		{
			FilmDim.y = FilmDim.x * ((real32)Image.Height/(real32)Image.Width);
		}
		else if(Image.Height > Image.Width)
		{
			FilmDim.x = FilmDim.y * ((real32)Image.Width/(real32)Image.Height);
		}
		
		int HalfImageWidth = Image.Width/2; 
		int HalfImageHeight = Image.Height/2; 
		
		uint32* Pixels = Image.Pixels;
		for(int y = -HalfImageHeight; y < HalfImageHeight; ++y)
		{
			real32 FilmY = (y*FilmDim.y/Image.Height);
			for(int x = -HalfImageWidth; x < HalfImageWidth; ++x)
			{
				real32 FilmX = (x*FilmDim.x/Image.Width);
				v3 Origin = V3(-2, 2, 0);
				v3 Dir = V3(Rotate(Mat4(), V3(0.5f, 1, 0), Radians(30)) * 
							V4(FilmX, FilmY, 1));
				
				v3 Color = RayTrace(Origin, Dir, &World, 1.0f, FLT_MAX, 3);
				// With gamma correction
				// v4 BMPValue = 
				// {
					// 255.0f*ExactLinearTosRGB(Color.x),
					// 255.0f*ExactLinearTosRGB(Color.y),
					// 255.0f*ExactLinearTosRGB(Color.z),
					// 255.0f
				// };
				// *Pixels++ = RGBAPack4x8(BMPValue);
				
				// Without gamma correction
				*Pixels++ = RGBAPack4x8(V4(Color * 255.0f, 255.0f));
			}
			
			if(y % 64 == 0)
			{
				printf("\rRayCasting... %d", (int)(100.0f * (real32)(y + HalfImageHeight) / (real32)Image.Height));
			}			
		}
		printf("\n");
		char* FileName = Args[1];
		WriteImage(FileName, &Image);
	}
	else
	{
		printf("Usage: %s width height file_name (To output an image)\n", Args[0]);
	}
	return 0;
}
